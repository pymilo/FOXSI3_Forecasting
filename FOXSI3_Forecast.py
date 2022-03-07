'''
Porpose:
      To download Farside images of the Sun to know ahead
      of time the solar activity for a particular day.
      Created for the FOXSI3 campaigne.
Action:
      Edit "MyPath" to set the location where you want to store
      images and reports.

Author: Milo Buitrago-Casas, SSL - UC Berkeley
Last update: August 13, 2018
'''

## Defining my Path:
MyPath = "/Users/Kamilobu/G-Drive/FOXSI/FOXSI-3 Flight 36.325/LaunchCampaign/FOXSI_Forecast"


## Packages for Sunpy processing : 
import sunpy.map
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from sunpy.physics.differential_rotation import diffrot_map, diff_rot
from astropy.coordinates import SkyCoord
import sunpy.map
import sunpy.coordinates
import sunpy.coordinates.wcs_utils
from sunpy.net import Fido, attrs as a
import urllib.request
## Packages for PDF Report:
from datetime import datetime, date, timedelta
from reportlab.platypus import SimpleDocTemplate
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import Image, Paragraph, Spacer, PageBreak
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY

## timing the program
Tinit = datetime.today()


''' Download Farside images of the Sun '''

urls = ['https://stereo-ssc.nascom.nasa.gov/where/where_is_stereo.gif',
        'https://stereo-ssc.nascom.nasa.gov/beacon/euvi_195_carrington_500x200.gif',
        'http://gong2.nso.edu/hot/Calibrated_Farside_Map.jpg',
        'http://farside.nso.edu/HMI/latest_cal_map.jpg',
        'https://stereo-ssc.nascom.nasa.gov/beacon/euvi_195_heliographic.gif',
        'http://jsoc.stanford.edu/data/farside/Recent/Composite_Map.png',
        'http://jsoc.stanford.edu/data/farside/Recent/AR_Map.png',
        'https://stereo.gsfc.nasa.gov/beacon/latest/ahead_euvi_195_latest.jpg',
        'http://jsoc.stanford.edu/data/farside/Recent/Movie.gif']
fnames = ['where_is_stereo.png',
          'beacon.png',
          'gongfarside.png',
          'hmifarside.png',
          'beacon_aia195.png',
          'Composite_Map.png',
          'AR_Map.png',
          'euvilatest.png',
          'movie.gif']

for ur,f in zip(urls,fnames):
  print('Downloading :'+f)
  urllib.request.urlretrieve(ur, MyPath+'/images/'+f)

''' Download latest Fits Files :  '''

files = ['fblos.fits','f0094.fits','f0171.fits','f0211.fits']

for file in files:
    print('Downloading :',file)
    url = 'http://sdowww.lmsal.com/sdomedia/SunInTime/'+datetime.utcnow().strftime("%Y/%m/%d")+'/'+file 
    urllib.request.urlretrieve(url, MyPath+'/fits/'+file)

''' Creating Differential Rotation Curve '''

latitudes = u.Quantity(np.arange(0, 90, 1), 'deg')
dt = 1 * u.day
rotation_rate = [diff_rot(dt, this_lat) / dt for this_lat in latitudes]
rotation_period = [360 * u.deg / this_rate for this_rate in rotation_rate]

fig = plt.figure(figsize=(10,7))
plt.plot(latitudes, [this_period.value for this_period in rotation_period])
plt.ylim(38, 24)
plt.ylabel('Rotation Period [{0}]'.format(rotation_period[0].unit))
plt.xlabel('Sin(Latitude)')
plt.title('Solar Differential Rotation Rate')
fig.savefig(MyPath+'/images/'+'diff_rot.png')


''' Defining FOXSI times: '''

TFOXSI = datetime(2018,9,7,12,0)
dt = ((TFOXSI - datetime.today()).total_seconds()/86400.) * u.day

''' Reading Fits files '''

latestHMIB = sunpy.map.Map(MyPath+'/fits/'+files[0])
latestAIA094 = sunpy.map.Map(MyPath+'/fits/'+files[1])
latestAIA171 = sunpy.map.Map(MyPath+'/fits/'+files[2])
latestAIA211 = sunpy.map.Map(MyPath+'/fits/'+files[3])
latestAIAFeXVIII = sunpy.map.Map(MyPath+'/fits/'+files[2])
latestAIAFeXVIII.data[:] = latestAIA211.data[:]/-120 + latestAIA171.data[:]/-450

## Rotate images:

rothmib = diffrot_map(latestHMIB,dt=dt)
rot094  = diffrot_map(latestAIA094,dt=dt)
rot171  = diffrot_map(latestAIA171,dt=dt)
rot211 = diffrot_map(latestAIA211,dt=dt)
## Temperate
TempMap = diffrot_map(latestAIA171,dt=dt)
TempMap.data[:] = rot211.data[:]/-120 + rot171.data[:]/-450

## HMI Magnetograms :

figHMIB = plt.figure(figsize=(15,7))
ax1HMIB = figHMIB.add_subplot(1, 2, 1)
latestHMIB.plot(vmin=-50,vmax=50,title='HMI Magnetogram Today, '+latestHMIB.date.strftime("%Y-%m-%d %H:%MUT"))
latestHMIB.draw_limb(color='red')
ax2HMIB = figHMIB.add_subplot(1, 2, 2)
rothmib.plot(vmin=-50,vmax=50,title='HMI Magnetogram rotated to '+TFOXSI.strftime("%Y-%m-%d %H:%M MDT"))
rothmib.draw_limb(color='red')
figHMIB.savefig(MyPath+'/images/HMIB.png')

## AIA 94 :

figAIA094 = plt.figure(figsize=(15,7))
ax1AIA094 = figAIA094.add_subplot(1, 2, 1)
latestAIA094.plot(vmin=0.0,vmax=1e2,title='AIA 094 Today, '+latestAIA094.date.strftime("%Y-%m-%d %H:%MUT"))
latestAIA094.draw_limb()
ax2AIA094 = figAIA094.add_subplot(1, 2, 2)
rot094.plot(vmin=0.0,vmax=1e2,title='AIA 094 rotated to '+TFOXSI.strftime("%Y-%m-%d %H:%M MDT"))
rot094.draw_limb()
figAIA094.savefig(MyPath+'/images/rot094.png')


## AIA 171 :

figAIA171 = plt.figure(figsize=(15,7))
ax1AIA171 = figAIA171.add_subplot(1, 2, 1)
latestAIA171.plot(title='AIA 171 Today, '+latestAIA171.date.strftime("%Y-%m-%d %H:%MUT"))
latestAIA171.draw_limb()
ax2AIA171 = figAIA171.add_subplot(1, 2, 2)
rot171.plot(title='AIA 171 rotated to '+TFOXSI.strftime("%Y-%m-%d %H:%M MDT"))
rot171.draw_limb(color='k')
figAIA171.savefig(MyPath+'/images/rot171.png')

## AIA Fe XVIII :

figAIAFeXVIII = plt.figure(figsize=(15,7))
ax1AIAFeXVIII = figAIAFeXVIII.add_subplot(1, 2, 1)
latestAIAFeXVIII.plot(vmin=-10,title='AIA FeXVIII Today, '+latestAIAFeXVIII.date.strftime("%Y-%m-%d %H:%MUT"))
latestAIAFeXVIII.draw_limb()
ax2AIAFeXVIII = figAIAFeXVIII.add_subplot(1, 2, 2)
TempMap.plot(vmin=-10,title='AIA Fe XVIII rotated to '+TFOXSI.strftime("%Y-%m-%d %H:%M MDT"))
TempMap.draw_limb(color='k')
figAIAFeXVIII.savefig(MyPath+'/images/rotFe18.png')


''' STEREO Downloading '''


'''
tstart = (datetime.today() - timedelta(days=6.307)).strftime("%Y-%m-%dT%H:%M:%S")
tend = (datetime.today() - timedelta(days=6.3)).strftime("%Y-%m-%dT%H:%M:%S")

stereo = (a.vso.Source('STEREO_A') &
          a.Instrument('EUVI') &
          #a.Time('2011-01-01', '2011-01-01T00:10:00'))
          a.Time(tstart,tend))
          
aia = (a.Instrument('AIA') &
       a.vso.Sample(24 * u.hour) &
       a.Time(tstart,tend))

wave = a.Wavelength(17 * u.nm, 18 * u.nm)
res = Fido.search(wave, aia | stereo)
print(res)

## Download
files = Fido.fetch(res)
print(files)

maps = {m.detector: m.submap(SkyCoord([-1100, 1100], [-1100, 1100],
                                      unit=u.arcsec, frame=m.coordinate_frame))
        for m in sunpy.map.Map(files)}

## Change coordinates : 

r = maps['AIA'].rsun_obs - 1 * u.arcsec  # remove the one arcsec so it's on disk.
# Adjust the following range if you only want to plot on STEREO_A
th = np.linspace(-180*u.deg, 0*u.deg)
x = r * np.sin(th)
y = r * np.cos(th)

coords = SkyCoord(x, y, frame=maps['AIA'].coordinate_frame)

r2 = maps['EUVI'].rsun_obs + 1 * u.arcsec  # remove the one arcsec so it's on disk.
# Adjust the following range if you only want to plot on STEREO_A
th2 = np.linspace(-180*u.deg, 0*u.deg)
x2 = r2 * np.sin(th2)
y2 = r2 * np.cos(th2)

coords2 = SkyCoord(x2, y2, frame=maps['EUVI'].coordinate_frame)


## Plot STEREO/EUVI and SDO/AIA : 

fig = plt.figure(figsize=(15, 7))
ax1 = fig.add_subplot(1, 2, 1, projection=maps['EUVI'])
maps['EUVI'].plot(axes=ax1)
ax1.plot_coord(coords, color='w')

ax2 = fig.add_subplot(1, 2, 2, projection=maps['AIA'])
maps['AIA'].plot(axes=ax2)
maps['AIA'].draw_limb()
fig.savefig(MyPath+'/images/'+'EUVI_AIA.png')
'''


''' ********************************************************* '''
'''                   Generating PDF Report                   '''

''' ********************************************************* '''



doc = SimpleDocTemplate(MyPath+"/FOXSI_Forecast_"+datetime.today().strftime("%Y%m%d_%H%M")+".pdf",
#doc = SimpleDocTemplate("test.pdf",
                        pagesize=landscape(letter),
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18)

Story=[]

styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
ptext = '<font size=18>%s</font>' % "Solar Forecast for the FOXSI-3 Launch"
Story.append(Paragraph(ptext, styles["Heading3"]))

Story.append(Spacer(1,20))

logo = MyPath+'/images/'+"FOXSI-3_small.png"
imlogo = Image(logo, 4*inch, 4*inch)
Story.append(imlogo)

ptext = "Date when this forecast was generated : "+datetime.today().strftime("%Y-%m-%d at %H:%M MDT")
Story.append(Paragraph(ptext, styles["Normal"]))

Story.append(Spacer(1,1))

ptext = '<font size=12>%s</font>' % "This forecast is produced for a final time of : "+TFOXSI.strftime('%b %d, %Y at %H:%Mpm')
Story.append(Paragraph(ptext, styles["Heading3"]))

Story.append(Spacer(1,20))

ptext = '<font size=10>%s</font>' % "Milo Buitrago-Casas"
Story.append(Paragraph(ptext, styles["Normal"]))       
ptext = '<font size=10>%s</font>' % "Space Sciences Laboratory"
Story.append(Paragraph(ptext, styles["Italic"]))
ptext = '<font size=10>%s</font>' % "UC Berkeley"
Story.append(Paragraph(ptext, styles["Italic"]))

#'''New Page'''
#Story.append(PageBreak())

#ptext = '<font size=12>%s</font>' % "Differential rotation curve"
#Story.append(Paragraph(ptext, styles["Heading1"]))
#diff_rot = MyPath+"/images/diff_rot.png"
#imdiff_rot = Image(diff_rot, 7*inch, 5*inch)
#Story.append(imdiff_rot)
#Story.append(Spacer(2,2))
#ptext = '<font size=12>%s</font>' % '"A COMPARISON OF DIFFERENTIAL ROTATION MEASUREMENTS", \
#Beck 1999 Solar Physics 191, 47–70. https://link.springer.com/content/pdf/10.1023%2FA%3A1005226402796.pdf'
#Story.append(Paragraph(ptext, styles["Normal"]))

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "HMI/SDO Magnetogram"
Story.append(Paragraph(ptext, styles["Heading1"]))
HMIB = MyPath+"/images/HMIB.png"
imhmib = Image(HMIB, 12*inch, 5.6*inch)
Story.append(imhmib)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "AIA/SDO 094"
Story.append(Paragraph(ptext, styles["Heading1"]))
aia094 = MyPath+"/images/rot094.png"
imaia094 = Image(aia094, 12*inch, 5.6*inch)
Story.append(imaia094)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "AIA/SDO 171"
Story.append(Paragraph(ptext, styles["Heading1"]))
aia171 = MyPath+"/images/rot171.png"
imaia171 = Image(aia171, 12*inch, 5.6*inch)
Story.append(imaia171)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "AIA/SDO Fe XVIII"
Story.append(Paragraph(ptext, styles["Heading1"]))
aiaFe18 = MyPath+"/images/rotFe18.png"
imaiaFe18 = Image(aiaFe18, 12*inch, 5.6*inch)
Story.append(imaiaFe18)

#'''New Page'''
#Story.append(PageBreak())

#ptext = '<font size=12>%s</font>' % "Where is STEREO-A right now?"
#Story.append(Paragraph(ptext, styles["Heading1"]))
#Story.append(Spacer(1,1))
#ptext = '<font size=12>%s</font>' % " STEREO has two separate telemetry\
#streams coming down from each spacecraft, the space weather beacon telemetry,\
#and the science recorder playback telemetry. The beacon telemetry contains \
#the most recent data and images, and is transmitted 24 hours per day. A volunteer\
#network of antenna stations around the world collect as much as possible of this \
#real-time data stream, and send it to the STEREO Science Center for processing. \
#However, because the beacon telemetry rate is very low, the images need to be \
#compressed by large factors, and are thus of much lower quality than the actual science data."
#Story.append(Paragraph(ptext, styles["Normal"]))
#Story.append(Spacer(5,5))
#ptext = '<font size=12>%s</font>' %"The science data collected by the STEREO spacecraft \
#are written to the on-board recorder, which is then read out and transmitted to the \
#ground during daily telemetry tracks using the NASA Deep Space Network. These data \
#are of much higher quality than the beacon data, but take several days to reach the \
#STEREO Science Center website. Thus, the most recent images on the STEREO Science \
#Center browse tool will always be beacon images. These temporary beacon images are replaced \
#with the full-quality versions as they become available, generally about 2-3 days later."
#whereS = MyPath+"/images/where_is_stereo.png"
#Story.append(Paragraph(ptext, styles["Normal"]))
#Story.append(Spacer(6,6))
#imwhereS = Image(whereS, 5*inch, 4*inch)
#Story.append(imwhereS)


'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "Latest EUVI 195 Beacon"
Story.append(Paragraph(ptext, styles["Heading1"]))
#Story.append(Spacer(1,1))
#ptext = '<font size=12>%s</font>' % "Shown here is the latest SECCHI \
#beacon image. The STEREO space weather beacon telemetry mode is a very \
#low rate, highly compressed data stream broadcast by the spacecraft 24 \
#hours per day. These data are used for space weather forecasting. Because \
#of the large compression factors used, these beacon images are of much lower \
#quality than the actual science data."
#Story.append(Paragraph(ptext, styles["Normal"]))
Story.append(Spacer(5,5))
euvilatest = MyPath+"/images/euvilatest.png"
imeuvilatest = Image(euvilatest, 6*inch, 6*inch)
Story.append(imeuvilatest)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "Beacon EUVI-195 & AIA/SDO-195"
Story.append(Paragraph(ptext, styles["Heading1"]))
beacon = MyPath+"/images/beacon_aia195.png"
imbeacon = Image(beacon, 10*inch, 6*inch)
Story.append(imbeacon)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "Latest EUVI-171 & AIA/SDO-171"
Story.append(Paragraph(ptext, styles["Heading1"]))
euviaia = MyPath+"/images/EUVI_AIA.png"
imeuviaia = Image(euviaia, 12*inch, 5.6*inch)
Story.append(imeuviaia)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "Composite & Active Regions Maps - Charlie Lindsey "
Story.append(Paragraph(ptext, styles["Heading1"]))
Composite_Map = "images/Composite_Map.png"
imComposite_Map = Image(Composite_Map, 7*inch, 3*inch)
Story.append(imComposite_Map)
AR_Map = "images/AR_Map.png"
imAR_Map = Image(AR_Map, 8*inch, 3*inch)
Story.append(imAR_Map)

'''New Page'''
Story.append(PageBreak())

ptext = '<font size=12>%s</font>' % "Beacon and GONG/HMI Farside "
Story.append(Paragraph(ptext, styles["Heading1"]))
beacon = MyPath+"/images/beacon.png"
imbeacon = Image(beacon, 5*inch, 2*inch)
Story.append(imbeacon)
gongfarside = MyPath+"/images/gongfarside.png"
imgongfarside = Image(gongfarside, 5*inch, 2*inch)
Story.append(imgongfarside)
hmifarside = MyPath+"/images/hmifarside.png"
imhmifarside = Image(hmifarside, 5*inch, 2*inch)
Story.append(imhmifarside)


euvilatest

doc.build(Story)

print('Listo el Pollo!!!')
print('Running time : ',(datetime.today() - Tinit).seconds/60.,' minutes.')


