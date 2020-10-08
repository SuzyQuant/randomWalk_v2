#!/usr/bin/env python

import matplotlib
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

data = np.loadtxt('quantitiesVStime2.txt')

time = data[:, 0] + 1
MSDth = time			#<x^2> analytical ; ALSO THE SECOND COLUMN
MSDnum = data[:, 1]		#<x^2> numerical / computer simulated
MSDerr = data[:, 3]		#error bar of <x^2>

Fsq1 = data[:, 4]		#F_s(q1,x(t))
Fsq1num = data[:, 10]		#F_s(q1,x(t)) analytical
Fsq1err = data[:, 7]	#error bar of F_s(q1,x(t))

Fsq2 = data[:, 5]		#F_s(q2,x(t))
Fsq2num = data[:, 11]		#F_s(q1,x(t)) analytical
Fsq2err = data[:, 8]	#error bar of F_s(q2,x(t))

Fsq3 = data[:, 6]		#F_s(q3,x(t))
Fsq3num = data[:, 12]		#F_s(q1,x(t)) analytical
Fsq3err = data[:, 9]	#error bar of F_s(q2,x(t))

with PdfPages('RESULTS20106.pdf') as pdf:

# ======== 1 log-log <x^2> VS time ==============	
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.title('log-log plot: <x^2> VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')

	plt.plot( time , MSDnum, 'r-*',label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 2 log-log <x^2> with error bars VS time ==============		
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.title('log-log plot w/ error bars: <x^2> VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.errorbar(time , MSDnum, MSDerr, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 3 log-log Fsq1 with error bars VS time ==============		
	plt.xlim([1,500])
	#plt.ylim([10**(-7),1.9])
	plt.title('log-log plot w/ error bars: Fsq1 VS time')
	plt.xlabel(r'$t$', fontsize=8)
	plt.ylabel(r'$log F_s(q1, t)$', fontsize=8)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.errorbar(time , Fsq1, Fsq1err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq1num, 'c-',label= r'Theory', linewidth=1)
	plt.tight_layout()

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 4 log-log Fsq2 with error bars VS time ==============		
	plt.xlim([1,500])
	plt.ylim([10**(-7),1.9])
	plt.title('log-log plot w/ error bars: Fsq2 VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$log F_s(q2, t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.errorbar(time , Fsq2, Fsq2err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq2num, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 5 log-log Fsq3 with error bars VS time ==============		
	plt.xlim([1,500])
	#plt.ylim([10**(-7),1.9])
	plt.title('log-log plot w/ error bars: Fsq3 VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$log F_s(q3, t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.errorbar(time , Fsq3, Fsq3err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq3num, 'c-',label= r'Theory', linewidth=1)	
	plt.tight_layout()

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== **** NO LOG PLOTS **** ==============	
# ======== 6 <x^2> VS time ==============	
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.title('<x^2> VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	
	plt.plot( time , MSDnum, 'r-*',label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 7 <x^2> with error bars VS time ==============		
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.title('w/ error bars: <x^2> VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	
	plt.errorbar(time , MSDnum, MSDerr, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 8 Fsq1 with error bars VS time ==============		
	plt.xlim([1,500])
	plt.ylim([-0.5,1.2])
	plt.title('w/ error bars: Fsq1 VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$F_s(q1, t)$', fontsize=10)
	
	plt.errorbar(time , Fsq1, Fsq1err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq1num, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 9 Fsq2 with error bars VS time ==============		
	plt.xlim([1,500])
	plt.ylim([-0.5,0.5])
	plt.title('w/ error bars: Fsq2 VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$F_s(q2, t)$', fontsize=10)
	
	plt.errorbar(time , Fsq2, Fsq2err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq2num, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======== 10 Fsq3 with error bars VS time ==============		
	plt.xlim([1,500])
	plt.ylim([-0.5,0.5])
	plt.title('w/ error bars: Fsq3 VS time')
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$F_s(q3, t)$', fontsize=10)
	
	plt.errorbar(time , Fsq3, Fsq3err, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , Fsq3num, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()
