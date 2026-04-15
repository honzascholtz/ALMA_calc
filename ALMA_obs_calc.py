'''
Based on script from Gareth Jones 
Modified by Jan Scholtz
'''

#Given SFR, redshift, velocity width, and SB/SFG classification, outputs expected flux density of CO lines.

import numpy as np
from decimal import Decimal
import sys
from math import *
from astropy.cosmology import Cosmology


class ALMA_calc:
	def __init__(self, gal):
		self.whichone = gal['NAME']
		self.sfr = gal['SFR']
		self.redsh = gal['ZRED']
		self.MSorSB = gal['MSorSB']
		self.alpha = gal['alpha']
		self.FWHM = gal['FWHM']
		self.numchan = gal['numchan']
		self.flipper = gal['flipper']
		self.ZP = gal['ZP']
		self.Cosmo = gal['Cosmo']
		self.gal= gal
		self.DL = self.GetDL( self.Cosmo)
		self.c=self.FWHM/(2*np.sqrt(2*np.log(2)))


	#Given SFR and galaxy class (MS or starburst), output expected MH2 - Sargent+14 e4: https://iopscience.iop.org/article/10.1088/0004-637X/793/1/19/pdf
	def SFR2MH2(self,):
		if self.MSorSB=='MS':
			alpha=9.22; beta=0.81 #MS
		else:
			alpha=8.05; beta=0.81 #SB
		logmh2=alpha+beta*np.log10(self.sfr)
		return 10**logmh2

	#Given MH2 and alpha_CO, returns L'CO(1-0)
	def MH2toLPRIMECO(self,):
		self.LPRIMECO = self.MH2/self.alpha
		print('alpha', self.alpha)
		return self.MH2/self.alpha

	#Narayanan+12
	def getLPCO10(self,ACO):
		temp1 = (ACO**-0.32) * (self.ZP**0.65) * self.MH2 / 10.7
		LP_N12 = temp1**(1/0.68)
		test_this = 10.7 * (LP_N12/ACO)**-0.32
		if test_this<6.3:
			N12a = 10.7 * ((LP_N12/ACO)**-0.32) / (self.ZP**0.65)
			print('N12 alpha:',N12a)
			return LP_N12	
		else:
			return (self.ZP**0.65) * self.MH2 / 6.3


	#Given L'CO(1-0), returns Sdv - Carilli & Walter (2013), p9: https://arxiv.org/pdf/1301.0371.pdf
	def LPRIMECOtoSDV(self,rj, nurest):
				   
		return ((1+self.redsh)**3)*((nurest/(1+self.redsh))**2)*(self.LPRIMECO10*rj)/((3.25E+7)*self.DL*self.DL)

	#1D Gaussian for spectrum fit
	def oneD_Gaussian(self,x,a,b,c):
		val1 = a * np.exp(-(x-b)**2 / (2*c*c))
		return val1

	#Given integral of Gaussian and the width, returns amplitude
	def getA(self,SDV,C):
		return SDV/(C*np.sqrt(2*np.pi))	

	#Integrates Gaussian over channels, returns expected flux [mJy]
	def getIC(self,A,C,numchan):
		temp=np.zeros(numchan)
		borders=np.linspace(-0.5,0.5,numchan+1)
		for CHAN in range(numchan):
			intrange=np.linspace(borders[CHAN]*C*(2*np.sqrt(2*np.log(2))),borders[CHAN+1]*C*(2*np.sqrt(2*np.log(2))),100)
			for i in intrange:
				temp[CHAN]+=self.oneD_Gaussian(i,A,0,C)
			temp[CHAN]=temp[CHAN]/len(intrange)
		return temp*1000.

	def SFRtoSDVCII(self,SFR, nucii, DL,LCII_Lsol=None ):
		if LCII_Lsol ==None:
			slope=1.18;offset=8.52 #high-z
			temp1=10**((np.log10(SFR)+offset)/slope)
			print('LCII based on SFR', np.log10(temp1))
			#slope=0.80;offset=5.73 #low-Z dwarf
		else:	
			temp1 = LCII_Lsol
			
		temp2=(1.04E-3)*(DL**2)*nucii
		return temp1/temp2

	def SFRtoSDVOIII(self,SFR, nuobs, DL, LOIII_Lsol=None):
		if LOIII_Lsol==None:
			LOIII_Lsol = 10**(7.4 + 0.97*np.log10(SFR)) #https://iopscience.iop.org/article/10.3847/1538-4357/ab94bd
			print('LOIII based on SFR', np.log10(LOIII_Lsol))
		#LOIII_Lsol = 10**(7.23 + 1.04*np.log10(SFR)) #Arata+21
		temp1 = (1.04E-3) * (DL**2) * nuobs
		return LOIII_Lsol/temp1

	#In Mpc
	def GetDL_gareth(self,):
		c = 299792.458
		H0=67.8; h=H0/100.
		az = 1.0/(1+1.0*self.zred)
		WK=0; WM=0.308; WR=4.165E-5/(h*h); WV=1-WM-WR-WK
		n=1000         # number of points in integrals
		DCMR = 0.0
		az = 1.0/(1+1.0*self.zred)
		a=1.0
		for i in range(n):
			a = az+(1-az)*(i+0.5)/n
			adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
			DCMR += 1./(a*adot)
		DCMR = (1.-az)*DCMR/n
		ratio = 1.00
		x = sqrt(abs(WK))*DCMR
		if x > 0.1:
			if WK > 0:
				ratio =  0.5*(exp(x)-exp(-x))/x 
			else:
				ratio = np.sin(x)/x
		else:
			y = x*x
			if WK < 0: y = -y
			ratio = 1. + y/6. + y*y/120.
		DCMT = ratio*DCMR
		DA = az*DCMT
		DL = DA/(az*az)
		DL_Mpc = (c/H0)*DL
		return DL_Mpc

	def GetDL(self, cosmo=None):
		if not cosmo:
			from astropy.cosmology import FlatLambdaCDM
			cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
		
		return cosmo.luminosity_distance(self.redsh).value


	def getSPWs(self,line_nu, flipper):
		nu1 = round(line_nu,3)
		nu2 = round(nu1 + 1.875,3) # - (5*0.03125)
		if line_nu>211 and line_nu<275:
			dSB = 10 + 4
		else:
			dSB = 8 + 4
		if flipper:
			nu3 = nu1 - dSB
			nu4 = nu2 - dSB
		else:
			nu3 = nu1 + dSB
			nu4 = nu2 + dSB
		return [nu1,nu2,nu3,nu4]


	def CII_calc(self, LCII_Lsol=None):
		self.DL = self.GetDL( self.Cosmo)
		c=self.FWHM/(2*np.sqrt(2*np.log(2)))

		self.nucii=1900.5369/(1.+self.redsh)
		print('[CII]158um')
		self.SDV=self.SFRtoSDVCII(self.sfr, self.nucii, self.DL, LCII_Lsol=LCII_Lsol)
		print('Sdv:',self.SDV,'Jy km/s')
		print('Rep Freq:',self.nucii,'GHz')
		print("Peak Amplitude: ",round(1000.*self.getA(self.SDV,c),2)/self.gal['PSF_n']*self.gal['mu']," mJy")
		print("Channel width: ",round(self.FWHM/self.numchan,1)," km/s")
		print([round(x,2) for x in self.getIC(self.getA(self.SDV,c),c,self.numchan)]," mJy")
		print('SPW centres:',self.getSPWs(self.nucii,self.flipper[0]))
		print("---")

	def OIII_calc(self, LOIII_Lsol=None):
	

		c=self.FWHM/(2*np.sqrt(2*np.log(2)))

		self.nuoiii=3393.006244/(1.+self.redsh)
		print('[OIII]88um')
		self.SDV=self.SFRtoSDVOIII(self.sfr, self.nuoiii, self.DL, LOIII_Lsol=LOIII_Lsol)
		print('Sdv:',self.SDV,'Jy km/s')
		print('Rep Freq:',self.nuoiii,'GHz')
		print("Peak Amplitude: ",round(1000.*self.getA(self.SDV,c),2)/self.gal['PSF_n']*self.gal['mu']," mJy")
		print("Channel width: ",round(self.FWHM/self.numchan,1)," km/s")
		print([round(x,2) for x in self.getIC(self.getA(self.SDV,c),c,self.numchan)]," mJy")
		print('SPW centres:',self.getSPWs(self.nuoiii,self.flipper[1]))
		print("-----------------------------------------------------------")


	def CO_calc(self, CO_list_calc ):
		
		co_list = []
		co_list.append({'J':1, 'nurest':115.27120180, 'rj1':1.0})
		co_list.append({'J':2, 'nurest':230.53800000, 'rj1':0.8})
		co_list.append({'J':3, 'nurest':345.79598990, 'rj1':0.7})
		co_list.append({'J':4, 'nurest':461.04076820, 'rj1':0.5})
		co_list.append({'J':5, 'nurest':576.26793050, 'rj1':0.4})
		co_list.append({'J':6, 'nurest':691.47307630, 'rj1':0.3})
		co_list.append({'J':7, 'nurest':806.65180600, 'rj1':0.2})
		co_list.append({'J':8, 'nurest':921.79970000, 'rj1':0.1})
		co_list.append({'J':9, 'nurest':1036.91239300, 'rj1':0.07})

		co_agn_list = []
		co_agn_list.append({'J':1, 'nurest':115.27120180, 'rj1':1.0})
		co_agn_list.append({'J':2, 'nurest':230.53800000, 'rj1':0.99})
		co_agn_list.append({'J':3, 'nurest':345.79598990, 'rj1':0.97})
		co_agn_list.append({'J':4, 'nurest':461.04076820, 'rj1':0.87})
		co_agn_list.append({'J':5, 'nurest':576.26793050, 'rj1':0.69})
		co_agn_list.append({'J':6, 'nurest':691.47307630, 'rj1':0.55})
		co_agn_list.append({'J':7, 'nurest':806.65180600, 'rj1':0.45})
		co_agn_list.append({'J':8, 'nurest':921.79970000, 'rj1':0.3})
		co_agn_list.append({'J':9, 'nurest':1036.91239300, 'rj1':0.2})

		if self.gal['CO_ladder'] == 'SF':
			co_list = co_list
		else:
			co_list = co_agn_list

		if 'MH2' not in self.gal.keys():
			self.MH2= self.SFR2MH2()
		else:
			self.MH2 = self.gal['MH2']
		N12 = self.MH2toLPRIMECO()#self.getLPCO10(np.pi*(1E+3)**2)
		print("MH2: ","{:.2E}".format(Decimal(str(self.MH2)))," Msol")
		print("FWHM: ","{:.2E}".format(Decimal(str(self.FWHM)))," km/s")
		self.LPRIMECO10=N12#M
		print('N12:'+"{:.2E}".format(Decimal(str(N12))))
		print("---")
			
		for i in range(len(co_list)):
			if co_list[i]['J'] in CO_list_calc:
				print("L'CO("+str(co_list[i]['J'])+"-"+str(co_list[i]['J']-1)+"): ","{:.2E}".format(Decimal(str(self.LPRIMECO10*co_list[i]['rj1']))),' K km/s pc^2')

				self.SDV=self.LPRIMECOtoSDV(co_list[i]['rj1'],co_list[i]['nurest'])

				print("Sdv: ",round(self.SDV*1E+3,2)," mJy km/s")
				print("Rep Freq:",round(co_list[i]['nurest']/(1+self.redsh),3)," GHz")
				print("Peak amplitude: ",round(1000.*self.getA(self.SDV,self.c)/self.gal['PSF_n']*self.gal['mu'],3)," mJy")
				print("Channel width: ",round(self.FWHM/self.numchan,1)," km/s")
				print([round(x,3) for x in self.getIC(self.getA(self.SDV,self.c),self.c,self.numchan)]," mJy")
				print('SPW centres:',self.getSPWs(co_list[i]['nurest']/(1+self.redsh),False))
				print("---")