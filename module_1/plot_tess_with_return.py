# syntax:   python3 plot_tess.py "starname" [--sector n --filter pcenter,pwidth --fake period,depth --zoom time1,time2 --flare nsig --deflc ]

from scipy.signal import medfilt
import requests
from bs4 import BeautifulSoup
import numpy as np
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle
import sys
import argparse
from astropy.stats import mad_std
from astroquery.mast import Catalogs, Tesscut
import argparse
from pathlib import Path
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.timeseries import BoxLeastSquares as BLS

def main(starName = "GJ 480", filterperiods = None, fake = None, 
         zoom = None, deflc = False, flare = None, sector = None, 
         prange = (0.2,10), detrend = False):
    # TIC star matching settings
    radSearch = 5 / 3600.
    nstar = 5

    # Lomb-Scargle periodogram settings
    #pmin = 0.2
    #pmax = 20.0
    prange= prange.split(",")
    pmin = float(prange[0])
    pmax = float(prange[1])
    nperiod = 10000

    ndetrend = 15 # detrending timescale

    if len(sys.argv) < 2:
        print ("Usage: check_tess.py starname (in double quotes)")
    else:
        starName = sys.argv[1]
        print('Star = ',starName)

        catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")

        print( catalogData[:nstar]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType'] )
        ticid = int(catalogData[0]['ID'])
        print('TIC counterpart = ',ticid)
        Ra = catalogData[0]['ra']
        Dec = catalogData[0]['dec']

        radNearby = 21.*10/3600.
        catalogData = Catalogs.query_object(starName, radius = radNearby, catalog = "TIC")
        nearby = (catalogData['Tmag'] < 15)
        #print( catalogData[:nstar]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType'] )
        nearbyStars = list( map( lambda x,y:[x,y], catalogData[nearby]['ra'], catalogData[nearby]['dec'] ) )
        nearbyNames = list( catalogData[nearby]['ID'] )

        coord = SkyCoord(Ra, Dec, unit="deg")
        sector_table = Tesscut.get_sectors(coordinates=coord)
        if not sector_table:
            print('WARNING: TESS has not observed this star')
        else:
            print(sector_table)
            sectors = sector_table['sector']
            sectors = list(sectors.astype(int))
            if sector == None:
                sector = sectors[0]
                print('You did not specify sector: using sector',sector)
            else:
                sector = int(sector)
                if sector not in sectors:
                    print('WARNING: Target not observed in this sector.  Using default') 
                    sector = sectors[0]
                
        if sector:
            sector = int(sector)
            print('Analyzing sector ',sector)
            ddir1 = "https://archive.stsci.edu/missions/tess/tid/"
            ddir2 = "s"+str(sector).zfill(4)+"/0000/"
            ddir3 = ""+str(int(ticid/1.0e8)).zfill(4)+"/"
            ddir4 = str(int((ticid-1.0e8*int(ticid/1.0e8))/1.0e4)).zfill(4)+"/"
            ddir5 = str(int(ticid - 1.0e4*int(ticid/1.0e4))).zfill(4)+"/"

            url = ddir1+ddir2+ddir3+ddir4+ddir5
            ext_lc = "s_lc.fits"
            ext_tp = "s_tp.fits"
            page = requests.get(url).text
            soup = BeautifulSoup(page,'html.parser')
            links1 = [a['href'] for a in soup.find_all('a') if a.get('href').endswith(ext_lc)]
            links2 = [a['href'] for a in soup.find_all('a') if a.get('href').endswith(ext_tp)]

            #r = requests.head(lcfile)
            #if ((r.status_code != 200) | (not links1) | (not links2)): 

            if ((not links1) | (not links2)): 
                print('No lightcurve file on MAST')
            else:
                lcfile =  url + links1[0] 
                tpfile =  url + links2[0]
                
                fig = plt.figure(figsize = (12,6))

                lchdu = fits.open(lcfile)
                lchdu.writeto(starName.replace(" ","_")+'_S'+str(sector)+'.fits',overwrite=True)
                fits.info(lcfile)

                lightcurvefile = starName.replace(" ","_")+'_S'+str(sector)+'_lc.csv'

                lcheader = lchdu[0].header
                lcdata = lchdu[1].data
                lchdu[1].columns
                aperture = lchdu[2].data

                ival = np.arange(9)
                ival = ival[::-1]
                for i in ival:
                    flag = aperture/2**i
                    flag = flag.astype(int)
                    aperture = aperture - 2**i*flag
                    if i == 1: sapap = flag
                    if i == 2: backap = flag

                    aperturefile = starName.replace(" ","_")+'_S'+str(sector)+'_sapap.txt'
                    backgroundfile = starName.replace(" ","_")+'_S'+str(sector)+'_backap.txt'

                if Path(backgroundfile).is_file():
                    backap = np.loadtxt(backgroundfile)
                else:
                    np.savetxt(backgroundfile,backap.astype(int),fmt='%1u')

                if Path(aperturefile).is_file():
                    sapap = np.loadtxt(aperturefile)
                else:
                    np.savetxt(aperturefile,sapap.astype(int),fmt='%1u')

                #sapflux = lcdata['SAP_FLUX']
                #pdcflux = lcdata['PDCSAP_FLUX']
                quality = lcdata['QUALITY']
                time = lcdata['TIME']

                tphdu = fits.open(tpfile)
                tpf_data = tphdu[1].data
                images = np.array(tpf_data['FLUX'])
                mean_image = np.nanmean(images,axis=0)

                nx, ny = mean_image.shape
                xx = np.outer(np.ones(nx),np.arange(ny))
                yy = np.outer(np.arange(nx),np.ones(ny))

                apertured_image = mean_image*sapap
                xc = np.mean(apertured_image*xx)/np.mean(apertured_image)
                yc = np.mean(apertured_image*yy)/np.mean(apertured_image)

                r = np.sqrt((xx-xc)**2 + (yy-yc)**2)
                r = np.ndarray.flatten(r)
                f = np.ndarray.flatten(mean_image)
                indices = np.argsort(r)
                fc = np.ndarray.cumsum(f[indices])

                fig.add_subplot(242)
                plt.plot(r[indices],fc/max(fc))
                plt.xlabel('aperture radius (pixels)',fontsize=9)
                plt.ylabel('fractional flux',fontsize=9)
                plt.title('curve of growth',fontsize=10)
                plt.ylim([0,1])

                wcs = WCS(tphdu[2].header)
                starLoc = wcs.all_world2pix([[Ra,Dec]],0)  #Second is origin
                nearbyLoc = wcs.all_world2pix(nearbyStars[1:], 0)

                fig.add_subplot(241, projection = wcs)
                plt.title(starName+' mean image',fontsize=10)

                saturation = 1.*(mean_image > 1.0e5)
                allones = mean_image*0 + 0.7
                colormap = plt.cm.Reds
                rgbval = colormap(allones)
                rgbval[:,:,3] = saturation

                plot_image = mean_image*(mean_image < 1.0e5)
                plt.imshow(np.sqrt(plot_image-np.min(plot_image)), origin = 'lower', cmap = plt.cm.viridis)
                plt.imshow(rgbval)
                plt.xlabel('RA', fontsize=9)
                plt.ylabel('Dec', fontsize=9)
                plt.grid(axis = 'both', color = 'white', ls = 'solid')
                plt.plot([xc],[yc],'o',c='black',markersize=3)
                
                x0 = starLoc[0,0]
                y0 = starLoc[0,1]

                plt.scatter(starLoc[0,0], starLoc[0,1], s = 30, color = 'red')
                #print(nearbyLoc)
                #nearbyLoc = np.array(nearbyLoc)
                #print(nearbyLoc.shape)
                #for xval,yval,name in zip(nearbyLoc[0:,0],nearbyLoc[0:,1],nearbyNames[1:]):
                #    if ((np.abs(xval-x0) < ny/2.) & (np.abs(yval-y0) < nx/2.)):
                        
                #        plt.plot(xval,yval,'o',c='white',markersize=3)
                #        plt.text(xval+0.4,yval+0.4,name,fontsize=5,color='white')
                        
                
                y = np.linspace(0,ny,ny*10)
                x = np.linspace(0,nx,nx*10)
                xx, yy = np.meshgrid(x[:-1],y[:-1])

                f = lambda x,y: backap[int(x),int(y)]
                g = np.vectorize(f)
                z = g(xx[:-1],yy[:-1])
                z = np.flip(z.T,axis=0)
                dx = -0.5
                dy = -0.5

                plt.contour(z[::-1], [0.5], colors='blue', extent=[dy,max(y)+dy,dx,max(x)+dx])

                f =  lambda x,y: sapap[int(x),int(y)]
                g = np.vectorize(f)
                z = g(xx[:-1],yy[:-1])
                z = np.flip(z.T,axis=0)
                plt.contour(z[::-1], [0.5], colors='orange', extent=[dy,max(y)+dx,dx,max(x)+dx])

                #if deflc == False:
                sapflux = np.einsum('ijk,jk',images,sapap)
                sapflux = sapflux - np.einsum('ijk,jk',images,backap)*np.sum(backap)/np.sum(sapap)
                          

                
              
                good = np.where((quality == 0) & (np.roll(quality,-1) == 0))[0]
                bad_bits = np.array([1,2,3,4,5,6,8,10,12])
                value = 0
                for v in bad_bits:
                    value = value + 2**(v-1)
                bad_data = np.bitwise_and(quality, value) >= 1 
                fluxcent_col = lcdata['MOM_CENTR1']
                fluxcent_row = lcdata['MOM_CENTR2']
                distance = ((fluxcent_col-np.nanmean(fluxcent_col))**2 + (fluxcent_row-np.nanmean(fluxcent_row))**2)**(0.5)
                mom_dump = np.bitwise_and(quality, 2**5) >= 1 

                t0 = min(time)
                yplot = sapflux[good]
                xplot = time[good] - t0
                yplot = yplot/np.nanmedian(yplot)



                fig.add_subplot(243)

                with fits.open(lcfile, mode="readonly") as hdulist:
                    tess_bjds = hdulist[1].data['TIME']
                    sap_fluxes = hdulist[1].data['SAP_FLUX']
                    pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
                    qual_flags = hdulist[1].data['QUALITY']

                medpdc = 1.*pdcsap_fluxes/np.nanmedian(1.*pdcsap_fluxes[np.isfinite(pdcsap_fluxes)])
                lower = 1. - 4*np.nanstd(medpdc)
                upper = 1. + 4*np.nanstd(medpdc)
                plt.vlines(time[mom_dump]-t0, 0.5, 1.5, colors = 'r', linestyle='dotted')
                plt.ylim([lower,upper])
            
                plt.plot(tess_bjds-t0, medpdc, '.',markersize=1)
                where_gt0 = np.where(qual_flags > 0)[0]
                plt.plot(tess_bjds[where_gt0]-t0, medpdc[where_gt0], 'ro',markersize=2.5,label = 'quality flag set')
                plt.xlabel("days",fontsize=9)
                plt.legend(loc = 'upper left')
                plt.title('normalized PDC lightcurve',fontsize=10)

                t = tess_bjds - min(tess_bjds)
                good = np.where((qual_flags == 0) & ~np.isnan(pdcsap_fluxes))[0]
                t = t[good]
                tbjd = tess_bjds[good]
                d = distance[good]
                ffix = pdcsap_fluxes[good]  
                ffix = ffix/np.nanmedian(ffix)

                    
                
                
                #ffix = sapflux[good]



                if flare != None:
                    flarecut = np.nanmean(ffix) + float(flare)*np.nanstd(ffix)
                    good = np.where(ffix < flarecut)
                    t = t[good]
                    d = d[good]
                    ffix = ffix[good]

                if fake != None:
                    fake = fake.split(",")
                    periodfake = float(fake[0])
                    depthfake = float(fake[1])
                    durfake = 12.5/24.*(periodfake/365.24)**0.333
                    phase = t/periodfake - np.random.random()
                    phase = phase - phase.astype(int)
                    phase = phase + 1.*(phase < 0.)
                    fake = depthfake*(np.abs(phase - 0.5) < 0.5*durfake)
                    ffix = ffix - fake

                period = pmin*np.exp(np.log(pmax/pmin)*np.arange(nperiod)/float(nperiod))
                freqs = 1./period
                ls = LombScargle(t,ffix)
                pgram = ls.power(freqs)
                siglevel = ls.false_alarm_level(0.001,method='baluev')
                peakls = period[np.argmax(pgram)]

                fig.add_subplot(244)
                plt.xlim(pmin,pmax)
                plt.xscale('log')
                plt.ylim(0,1.2*max(pgram))
                plt.plot(period,pgram,'k',label='Lomb-Scargle')

                plt.title('Lomb-Scargle periodogram',fontsize=10)
                plt.plot([pmin,pmax],[siglevel,siglevel],'g')
                plt.xlabel('period (days)',fontsize=9)
                plt.plot([peakls],[1.05*max(pgram)],'ro',markersize=3)
                plt.text(1.3*pmin,1.07*max(pgram),'peak = '+str(round(peakls,3))+' days')
                plt.vlines(2.5,0,2*max(pgram),color='r',linestyle='dotted')
                plt.vlines(6.85,0,2*max(pgram),color='m',linestyle='dotted')
                plt.text(2.5,0.8*max(pgram),'mom dump',fontsize=7,rotation=90)
                plt.text(6.85,0.8*max(pgram),'half orbit',fontsize=7,rotation=90)
                print('peak periodic single from L-S analysis = ',peakls, ' days')

                peakls = 0.381076927  #### DELETE IF NOT K2-22 ####
                tc = 2456811.1207
                phase = (tbjd-tc)/peakls
                
                phase = phase - phase.astype(int)
                phase = phase - 1*(phase >= 0.5)
                phase = phase + 1.*(phase < -0.5)
                indices = np.argsort(phase)
                fphase = ffix[indices]
                ph = phase[indices]
                tph = t[indices]


                
                #f3 = np.concatenate((fphase,fphase,fphase))
                #xfit = np.arange(len(fphase))/(1.*len(fphase))
                nseg = 20
                fmed = ph*0
                nmc = 10000
                flow = 0*ph
                fhigh = 0*ph
                duration = 46./60/24. # duration of K2-22 transit
                durphase = duration/peakls
                
                in_transit = (np.abs(ph) < durphase/2.)
                out_transit = (np.abs(ph) > durphase/2.)
                fmedin = np.nanmedian(fphase[in_transit])
                fmedout = np.nanmedian(fphase[out_transit])
                fvals = []
                fin = fphase[in_transit]
                for imc in np.arange(nmc):
                    indices = np.random.randint(0,len(fin),len(fin))
                    fvals.append(np.nanmedian(fin[indices]))
                fvals = np.sort(np.array(fvals))
                flimit = fvals[int(0.01*nmc)]
                fupper = fvals[int(0.99*nmc)]
                print('measured and 99% lower limit on % transit depth =', 100*(fmedout-fmedin),100*(fmedout-flimit),100*(fmedout-fupper))
                
                nmc = 1000
                for iseg in np.arange(nseg):
                    which = ((ph >= (iseg/(1.*nseg)-0.5)) & (ph < ((iseg+1)/(1.*nseg))-0.5))
                    xfit = ph[which]
                    pfit = fphase[which]
                    #fpoly = np.polyfit(xfit,fphase[which],4)
                    #fmed[which] = fpoly[0]*xfit**4 + fpoly[1]*xfit**3 + fpoly[2]*xfit**2 + fpoly[3]*xfit + fpoly[4]
                    #fpoly = np.polyfit(xfit,fphase[which],2)
                    #fmed[which] = fpoly[0]*xfit**2 + fpoly[1]*xfit + fpoly[2]
                    fmed[which] = np.nanmedian(pfit)
                    fvals = []
                    for imc in np.arange(nmc):
                        indices = np.random.randint(0,len(pfit),len(pfit))
                        fvals.append(np.nanmedian(pfit[indices]))
                    
                    fvals = np.sort(np.array(fvals))
                    flow[which] = fvals[int(0.16*nmc)]
                    fhigh[which] = fvals[int(0.84*nmc)]
                    #print(fvals[int(0.16*nmc)],fvals[int(0.84*nmc)])
                        
                    
                #fmed = medfilt(f3,301)
                #fmed = fmed[len(fphase):2*len(fphase)]
                fh = open(lightcurvefile,'w+')
                for tval,phaseval,fph,fmedval,flowval,fhighval in zip(t,ph,fphase,fmed,flow,fhigh):
                    printline = str(tval)+','+str(phaseval)+','+str(fph)+','+str(fmedval)+','+str(flowval)+','+str(fhighval)+'\n'
                    fh.write(printline)
                
                fig.add_subplot(245)
                plt.plot(phase,ffix,'.',markersize=0.5)
                plt.plot(ph,fmed,'r')
                plt.plot(ph,flow,'r',alpha=0.5)
                plt.plot(ph,fhigh,'r',alpha=0.5)
                fcorr = fmed[np.argsort(tph)]
                if detrend:
                    fdetrend = ffix/fcorr
                else:
                    print('warning: not detrending the light curve!')
                    fdetrend = ffix*1.
                plt.ylim(0.99,1.01)
                plt.title('phased to to peak period',fontsize=10)
                plt.xlabel('phase',fontsize=9)

                durmin = 0.5*12.5/24.*(pmin/365.24)**0.333
                durmax = 2.0*12.5/24.*(pmax/365.24)**0.333
                durmax = durmax*(durmax < pmin) + 0.99*pmin*(durmax > pmin)
                durations = np.linspace(durmin,durmax,20)
                model = BLS(t,fdetrend)
                results = model.autopower(durations,frequency_factor=5.0, minimum_period=pmin,maximum_period=pmax)
                
                powerspec = results.power
                periodval = results.period

                #whichperiods = ((periodval > pmin) & (periodval < pmax))
                #powerspec = powerspec[whichperiods]
                #periodval = periodval[whichperiods]
                

                powerfilter = 1. + 0*powerspec
                if filterperiods != None:
                    filterperiods = filterperiods.split(",")
                    pfiltcenter = float(filterperiods[0])
                    pfiltwidth = float(filterperiods[1])
                    nfilter = int(filterperiods[2])
                    print('Filtering '+str(pfiltcenter)+' plus harmonics')
                    for n in range(1,nfilter+2):
                        powerfilter = powerfilter*(np.abs(periodval - pfiltcenter*n) > pfiltwidth/2.*n)
                        powerfilter = powerfilter*(np.abs(periodval - pfiltcenter/n) > pfiltwidth/2./n)
                    powerspec = powerspec*powerfilter

                fig.add_subplot(246)
                plt.plot(periodval, 1000*powerspec, "k", lw=0.5,label='BLS')
                index = np.argmax(powerspec)
                peakperiod = periodval[index]
                period = peakperiod
    
                plt.axvline(peakperiod, alpha=0.4, lw=3, color='green')
                for n in range(2, 10):
                    plt.axvline(n*peakperiod, alpha=0.6, lw=1, linestyle="dashed", color='green')
                    plt.axvline(peakperiod / n, alpha=0.6, lw=1, linestyle="dashed", color='green')
                ymax = 1300*max(powerspec)
                plt.xlim([min(periodval),max(periodval)])
                plt.vlines(2.5,0,ymax,color='r',linestyle='dotted')
                plt.vlines(6.85,0,ymax,color='m',linestyle='dotted')
                plt.vlines(13.7,0,ymax,color='m',linestyle='dotted')
                plt.ylim([0,ymax])
                plt.xlabel("period (days)",fontsize=9)
                plt.xscale('log')
                plt.ticklabel_format(axis='y',style='sci')
                
                if filterperiods != None:
                    for n in range(1,nfilter+2):
                        pfilt1 = pfiltcenter*n + pfiltwidth*n/2.
                        pfilt2 =  pfiltcenter*n - pfiltwidth*n/2.
                        plt.fill([pfilt1,pfilt2,pfilt2,pfilt1],[0,0,ymax,ymax],color='grey',alpha=0.5)  
                        pfilt1 = pfiltcenter/n + pfiltwidth/2./n
                        pfilt2 =  pfiltcenter/n - pfiltwidth*n/2./n
                        plt.fill([pfilt1,pfilt2,pfilt2,pfilt1],[0,0,ymax,ymax],color='grey',alpha=0.5)   

                plt.ylabel("power x 1000")
                plt.title('box least-squares search',fontsize=10)
                plt.text(1.2*min(periodval),1300*0.91*max(powerspec),'peak = '+str(round(peakperiod,3))+' days')

                tc0 = results.transit_time[index]
                duration = results.duration[index]
                depth = results.depth[index]

                modelstats = model.compute_stats(period, duration, tc0)
        
                print('Transit fit parameters = ', period,tc0,duration,depth)
                depthsig = 1/modelstats['depth'][1]
                print('Significance of BLS signal = ',depthsig)

                fig.add_subplot(247)
                x = (t - tc0 + 0.5*period) % period - 0.5*period
                m = np.abs(x) < 0.5 
                plt.plot(x[m]*24, 1000*(fdetrend[m]-1), ".k", ms=3)

                x = np.linspace(-12, 12, 1000) #* u.day
                f = model.model(x + tc0, period, duration, tc0)
                plt.plot(x*24, 1000*(f-1),c='r',linewidth=2)
                plt.xlim(-12, 12)
                plt.ylim(-5*depth*1000,3*depth*1000)
                plt.xlabel("time since transit (hr)",fontsize=9)
                plt.ylabel("differential signal (ppt)",fontsize=9)
                plt.title('transit lightcuve fit',fontsize=10)
                plt.text(-11,2400*depth,'depth = '+str(round(depth*1000,3))+' (ppt)')


                t1 = tc0 - 13
                t2 = tc0 + 13
                if zoom != None:
                    zoom = zoom.split(',')
                    t1 = float(zoom[0])
                    t2 = float(zoom[1])
                which = ((t > t1) & (t < t2) & (fdetrend > 0.))
            
                fig.add_subplot(4,4,12)
                #plt.ylim([1+0.8*(min(fdetrend[which])-1),1+1.2*(max(fdetrend[which])-1)])
                plt.ylim(min(fdetrend[which]),max(fdetrend[which]))
                plt.xlim([t1,t2])
                plt.vlines([tc0], 0.8*min(fdetrend[which]), 1.2*max(fdetrend[which]), colors = 'green', linestyle='dotted')
                plt.plot(t[which],fdetrend[which],'.',markersize=1)
                plt.vlines(time[mom_dump]-t0, 0.8*min(fdetrend[which]), 1.2*max(fdetrend[which]), colors = 'r', linestyle='dotted')
                fig.add_subplot(4,4,16)
                plt.xlim([t1,t2])
                plt.plot(t[which],d[which],'.',markersize=1)
                plt.ylim(0,max(d[which]))
                plt.vlines(time[mom_dump]-t0, 0, max(d[which]), colors = 'r', linestyle='dotted')
                plt.xlabel('days',fontsize=9)
                plt.ylabel('offset (pixels)',fontsize=9)

                plt.subplots_adjust(bottom=0.1, right=0.97, top=0.95, wspace=0.3, hspace=0.3, left=0.05)
                plotname = starName.replace(" ","_")+'_S'+str(sector)+'.png'
                plt.savefig(plotname)
                plt.show(block=True)

                return period, depth 

if __name__ == "__main__": 
    main()