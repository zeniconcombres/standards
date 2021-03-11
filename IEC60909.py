import math
from numpy import log as ln
#####IEC60909 Electric cables – Calculation of the current rating############

def Z1 (R_con,f,con_OD,s1,overall_OD,laying): ###Without metallic sheath
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	## Rc, conductor AC per km unitresistance
	omega=2*math.pi*f ###Angular frequence
	mu0=(4*math.pi)/(10**4) ### absolute permeability of vacuum
	A=omega*mu0/(2*math.pi)
	### laying format True in flat, False in trefoil
	if laying==True:
		d=s1*(2**(1/3))
	elif laying==False:
		d=overall_OD*1.06
	B=1/4+ln(d/(con_OD/2))
	X1=A*B
	Z1=complex(R_con,X1)
	return Z1

def Z1_s(Z1,f,Rs,overall_OD,s1,sheath_D,laying): ###With metallic sheath
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	### sheath_D Sheath mean diameter
	if laying==True:
		d=s1*(2**(1/3))
	elif laying==False:
		d=overall_OD*1.06
	omega=2*math.pi*f
	mu0=(4*math.pi)/(10**4)
	A=omega*mu0/(2*math.pi)
	B=ln(d/(sheath_D/2))
	C=A*B
	D=C**2
	E=complex(Rs,C)
	Z1_s=Z1+D/E
	return Z1_s

def Earth_p(R_soil,f):
	"""Requires the DC resistance input parameter to be in units Ohm.m
	for this to have the right scale of effect."""
	mu0=(4*math.pi)/(10**4)
	omega=2*math.pi*f
	A=(omega*mu0/R_soil)**0.5
	delta=1.851/A
	return delta


def Z0 (R_con,f,Earth_p,con_OD,overall_OD,s1, laying):
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	if laying==True:
		d=s1*(2**(1/3))
	elif laying==False:
		d=overall_OD*1.06
	mu0=(4*math.pi)/(10**4)
	omega=2*math.pi*f
	A=R_con+3*omega*mu0/8
	B=omega*mu0/(2*math.pi)
	C=((con_OD/2)*(d**2))**(1/3)
	D=1/4+3*ln(Earth_p*1000/C)
	E=B*D
	Z0=complex(A,E)
	return Z0

def Z0_s(Z0,f,Earth_p,sheath_D,overall_OD, Rs,s1, laying):
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	if laying==True:
		d=s1*(2**(1/3))
	elif laying==False:
		d=overall_OD*1.06
	mu0=(4*math.pi)/(10**4)
	omega=2*math.pi*f
	A=3*omega*mu0/8
	B=3*omega*mu0/(2*math.pi)
	C=((sheath_D/2)*(d**2))**(1/3)
	D=ln(Earth_p*1000/C)
	E=(complex(A,B*D))**2
	F=complex(Rs+A,B*D)
	Z0_s=Z0-E/F
	return Z0_s

def r3(Rs,s1,f,sheath_D,Earth_p,overall_OD,laying):
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	if laying==True:
		d=s1*(2**(1/3))
	elif laying==False:
		d=overall_OD*1.06
	mu0=(4*math.pi)/(10**4)
	omega=2*math.pi*f
	A=3*omega*mu0/8
	B=3*omega*mu0/(2*math.pi)
	C=ln(Earth_p*1000/(((sheath_D/2)*d**2)**(1/3)))
	D=complex((Rs+A),(B*C))
	r3=Rs/D
	return r3