#!/home/shanksj/misc/anaconda/bin/python

#import subprocess
import ctypes as C
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import re
import scipy.optimize as sciop
from matplotlib.widgets import Slider
import numdifftools as nd

flip_phase = 1

bmadlib = C.CDLL("/home/ehrlichm/openmp_bmad/production/lib/libphase_advance_4py.so")

bmadlib.parse_in_file.argtypes = [C.c_char_p, C.c_int, C.POINTER(C.c_int), C.POINTER(C.c_int), C.POINTER(C.c_char), \
                                  C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_char), \
																	C.POINTER(C.c_double), C.POINTER(C.c_double)]
bmadlib.init.argtypes = []
bmadlib.phase_vertical_calc.argtypes = [C.POINTER(C.c_double), C.c_double, C.c_double, C.POINTER(C.c_double)]

# Call external Fortran subroutine to parse .in file
in_file = b'phase_advance.in'
in_file_len = C.c_int(len(in_file))
c_nbpms = C.c_int(-1)
c_nparams = C.c_int(-1)
c_low_cur = C.c_double(-1.0)
c_high_cur = C.c_double(-1.0)
c_param_names = C.create_string_buffer(20000)
c_datafile = C.create_string_buffer(200)
upper_bounds = np.empty(1001, dtype=C.c_double) #MAX_PARAMS=10
lower_bounds = np.empty(1001, dtype=C.c_double) #MAX_PARAMS=10
initial_params = np.empty(1001, dtype=C.c_double) #MAX_PARAMS=10
bmadlib.parse_in_file(in_file,in_file_len,c_nbpms,c_nparams,c_param_names, upper_bounds.ctypes.data_as(C.POINTER(C.c_double)), \
                      lower_bounds.ctypes.data_as(C.POINTER(C.c_double)),initial_params.ctypes.data_as(C.POINTER(C.c_double)), \
											c_datafile, c_low_cur, c_high_cur)
param_names = c_param_names.value.decode("utf-8").split(',')[:-1]
datafile = c_datafile.value.decode("utf-8").strip()
nbpms = c_nbpms.value
nparams = c_nparams.value
low_cur = c_low_cur.value
high_cur = c_high_cur.value
upper_bounds = upper_bounds[0:nparams]
lower_bounds = lower_bounds[0:nparams]
bounds = np.column_stack((lower_bounds,upper_bounds))
initial_params = initial_params[0:nparams]


# Call external Fortran subroutine to parse lattice file and locate BPMs and parameters.
# Also initialize global module data.
bmadlib.init()

print("Number of BPMs: {}".format(nbpms))
print("Number of parameters: {}".format(nparams))

def load_data():
	nbpm = 0
	data = []
	with open(datafile,'r') as f:
		for l in f:
			if not l.startswith('#') and len(l.split()) > 0:
				data.append(l.split())
				nbpm += 1

	bpmNames = []
	phi_x = np.zeros(nbpm)
	phi_x_error = np.zeros(nbpm)
	phi_y = np.zeros(nbpm)
	phi_y_error = np.zeros(nbpm)

	i=0
	for l in data:
		bpmNames.append(l[0])
		phi_x[i] = float(l[1])
		phi_x_error[i] = float(l[2])
		phi_y[i] = float(l[3])
		phi_y_error[i] = float(l[4])
		i+=1

	return bpmNames,phi_x,phi_x_error,phi_y,phi_y_error

bpmNames, data_phi_x, data_phi_x_error, data_phi_y, data_phi_y_error = load_data()

print("FOO: ", bpmNames)
bpmNames_only = [x.split(":")[1] for x in bpmNames]
print("FOO C")

best_merit = 99999.9
best_p = np.zeros(nparams)

def evaluate(p,low_I,high_I):
	global best_merit, best_p

	delta_phi = np.empty(nbpms, dtype=C.c_double)
	bmadlib.phase_vertical_calc(p.ctypes.data_as(C.POINTER(C.c_double)),C.c_double(low_I),C.c_double(high_I),delta_phi.ctypes.data_as(C.POINTER(C.c_double)))

	#if len(bpmNames) != len(bpmNames_sim):
	#	print("Error: data file and simulation results have different number of BPMs!")
	#	sys.exit()		

	delta_phi = flip_phase*delta_phi

	residual = (delta_phi-data_phi_y)/data_phi_y_error
	merit = np.linalg.norm(residual,ord=2)
	rchi2 = merit**2/(nbpms-1)
	if(merit < best_merit):
		best_merit = merit
		best_p = p

	print('parameters: '+''.join('{:13.5e} '.format(x) for x in p)+'merit: {:13.5e} best: {:13.5e} rchi2: {:13.5e}'.format(merit, best_merit, rchi2))

	return merit,delta_phi

def evaluate_wrapper(p):
	merit,_ = evaluate(p,low_cur,high_cur)
	return merit

merit,sim_phi_y = evaluate(initial_params,low_cur,high_cur)
best_merit = merit

def parameter_uncertainty(params):
	step=np.abs(params)*0.001
	step[step<1.0e-6]=1.0e-6
	H = nd.Hessian(evaluate_wrapper,step=step)(params)
	cov=np.linalg.inv(H)
	# for i in cov:
	# 	print(''.join('{:13.5e} '.format(x) for x in i))
	return [np.sqrt(cov[i,i]) for i in range(len(cov))]

def drad_dI(p,punc):
	dradm = np.zeros_like(p)
	drad0 = np.zeros_like(p)
	dradp = np.zeros_like(p)
	for i in range(len(p)):
			pwrk = np.zeros_like(p)
			pwrk[i] = p[i] - punc[i]
			merit,sim_dphi = evaluate(pwrk,0.0,1.0)
			dradm[i] = sim_dphi[-1]

			pwrk = np.zeros_like(p)
			pwrk[i] = p[i]
			merit,sim_dphi = evaluate(pwrk,0.0,1.0)
			drad0[i] = sim_dphi[-1]

			pwrk = np.zeros_like(p)
			pwrk[i] = p[i] + punc[i]
			merit,sim_dphi = evaluate(pwrk,0.0,1.0)
			dradp[i] = sim_dphi[-1]
	
	return dradm, drad0, dradp

print("Calculate uncertainty of initial parameters? (y/[n]): ",end='')
if input() == 'y':
	punc = parameter_uncertainty(initial_params)

	dradm, drad0, dradp = drad_dI(initial_params,punc)
	dHzdim = dradm / (2.0*np.pi) / (2.56e-6)
	dHzdi0 = drad0 / (2.0*np.pi) / (2.56e-6)
	dHzdip = dradp / (2.0*np.pi) / (2.56e-6)
else:
	punc = np.zeros(len(initial_params))
	dHzdim = np.zeros(len(initial_params))
	dHzdi0 = np.zeros(len(initial_params))
	dHzdip = np.zeros(len(initial_params))

print("k1l and uncertainty")
for i in range(len(initial_params)):
	print('{0:20s}: {1:13.5e} +- {2:13.5e}'.format(param_names[i],initial_params[i],punc[i]))

print("Hz/mA and uncertainty range")
for i in range(len(initial_params)):
	print('{0:20s}: {1:13.5f}   {2:13.5f}   {3:13.5f}'.format(param_names[i],dHzdi0[i], dHzdi0[i]-dHzdim[i], dHzdi0[i]-dHzdip[i]))

x=list(range(len(data_phi_y)))

short_names = ["_".join(x.split("_")[0:-1]) for x in param_names]
print(short_names)

def sig_plot(p,param_names):
	matplotlib.rcParams.update({'font.size': 12})
	fig2, ax = plt.subplots(7,sharex=True,figsize=(20,12), gridspec_kw = {'height_ratios':[2,1,1,1,1,1,1]})
	[x.xaxis.grid(True) for x in ax]
	_,sim_phi_y = evaluate(p,low_cur,high_cur)
	ax[0].xaxis.tick_top()
	ax[0].xaxis.set_tick_params(labeltop=True)
	ax[0].set_xticks(x)
	ax[0].set_xticklabels(bpmNames_only, rotation='vertical')
	ax[0].set_xlabel("BPM Name")
	ax[0].plot(x,sim_phi_y*1000,label='combined effect of all simulated impedance sources')
	ax[0].errorbar(x,data_phi_y*1000,yerr=data_phi_y_error*1000,label='data',fmt='--o',capsize=5)
	ax[0].set_xlim(x[0]-1,x[-1]+1)
	ax[0].legend()
	for ix in range(len(p)-1): 
		p_ix = np.zeros_like(p)
		p_ix[ix] = p[ix]
		_,sim_phi_y = evaluate(p_ix,low_cur,high_cur)
		ax[ix+1].plot(x,sim_phi_y*1000,label='sim')
		ax[ix+1].set_xlim(x[0]-1,x[-1]+1)
	ax[-1].set_xticks(x)
	ax[-1].set_xticklabels(bpmNames_only, rotation='vertical')
	ax[-1].set_xlabel("BPM Name")
	ax[1].text(0.95,0.75,'scraper at 43W',bbox={'pad':5,'facecolor':'white'},transform=ax[1].transAxes,horizontalalignment='right')
	ax[2].text(0.95,0.75,'collimator at 43E',bbox={'pad':5,'facecolor':'white'},transform=ax[2].transAxes,horizontalalignment='right')
	ax[3].text(0.95,0.75,'undulator between 7W and 8W',bbox={'pad':5,'facecolor':'white'},transform=ax[3].transAxes,horizontalalignment='right')
	ax[4].text(0.95,0.75,'separators at 9W, 44W, 44E, and 8E',bbox={'pad':5,'facecolor':'white'},transform=ax[4].transAxes,horizontalalignment='right')
	ax[5].text(0.95,0.75,'RF cavities between 8E and 9E, and 8W and 9W',bbox={'pad':5,'facecolor':'white'},transform=ax[5].transAxes,horizontalalignment='right')
	ax[6].text(0.95,0.75,'RW & Bulk represented by 100 equally spaced moments of equal strength',bbox={'pad':5,'facecolor':'white'},transform=ax[6].transAxes,horizontalalignment='right')
	[a.set_ylabel(n) for a,n in zip(ax,['all sources', 'scraper','collimator','undulator','separators','RF cavities','RW & Bulk'])]
	fig2.text(0.06,0.5, r'$\phi_{high\ I}-\phi_{low\ I} \left(mrad\right)$', ha='center', va='center', rotation='vertical', fontsize='22')
	fig2.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in fig2.axes[1:-1]], visible=False)
	fig2.savefig('signatures.eps', format='eps') #, dpi=1000)

def pub_plot(p):
	matplotlib.rcParams.update({'font.size': 12})
	fig1 = plt.figure(figsize=(19,5))
	ax1 = fig1.add_axes([0.1,0.15,0.85,0.75,])
	ax1.xaxis.grid(True)
	_,sim_phi_y = evaluate(p,low_cur,high_cur)
	ax1.plot(x,sim_phi_y*1000,label='sim')
	ax1.errorbar(x,data_phi_y*1000,yerr=data_phi_y_error*1000,label='data',fmt='--o',capsize=5)
	ax1.set_xticks(x)
	ax1.set_xticklabels(bpmNames_only, rotation='vertical')
	ax1.set_xlabel("BPM Name")
	ax1.set_ylabel(r'$\phi_{high\ current}-\phi_{low\ current} \left(mrad\right)$', fontsize='18')
	ax1.set_ylim(-2.5,2.5)
	ax1.legend()
	fig1.savefig('fitter.eps', format='eps') #, dpi=1000)

slider_scale = 1.0e6
def interactive_plot(p):
	matplotlib.rcParams.update({'font.size': 10})
	fig = plt.figure(figsize=(19,10))
	ax = fig.add_subplot(111)
	ax.xaxis.grid(True)
	fig.subplots_adjust(left=0.10,bottom=0.40)
	_,sim_phi_y = evaluate(p,low_cur,high_cur)
	[sim_line] = ax.plot(x,sim_phi_y,label='sim')
	ax.errorbar(x,data_phi_y,yerr=data_phi_y_error,label='data',fmt='--o',capsize=5)
	[res_line] = ax.plot(x,data_phi_y-sim_phi_y,label='residual')
	plt.xticks(x, bpmNames, rotation='vertical')
	plt.legend()

	axis_color='lightgoldenrodyellow'
	init_slider=[]
	slider_height = 0.28/nparams
	for i in range(nparams):
		init_slider_ax = fig.add_axes([0.10, 0.27-slider_height*i, 0.75, slider_height], facecolor=axis_color)
		if param_names[i][0:2] == '__':
			init_slider.append(Slider(init_slider_ax, param_names[i], slider_scale*bounds[i][0], slider_scale*bounds[i][1], valinit=slider_scale*p[i], valfmt="%9.3f mrad"))
		else:
			init_slider.append(Slider(init_slider_ax, param_names[i], slider_scale*bounds[i][0], slider_scale*bounds[i][1], \
			                   valinit=slider_scale*p[i], valfmt="%9.3f x "+r"$10^{-6}$ k1l / mA"))

	def sliders_on_changed_init(val):
		p = np.zeros(nparams)
		for i in range(nparams):
			p[i] = init_slider[i].val/slider_scale
		merit,sim_phi_y = evaluate(p,low_cur,high_cur)
		sim_line.set_ydata(sim_phi_y)
		res_line.set_ydata(data_phi_y-sim_phi_y)
		fig.canvas.draw_idle()

	for i in range(nparams):
		init_slider[i].on_changed(sliders_on_changed_init)

	fig.show()

pub_plot(initial_params)
sig_plot(initial_params,param_names)
interactive_plot(initial_params)

print("Continue to optimization? (y/[n]): ",end='')
if input() is not 'y':
	sys.exit()

#res = sciop.minimize(evaluate_wrapper, initial_params, bounds=bounds)
#res = sciop.minimize(evaluate_wrapper, initial_params)
res = sciop.differential_evolution(evaluate_wrapper,bounds,disp=True)
#res = sciop.basinhopping(evaluate_wrapper,initial_params)

xopt = res.x

print()
print("----------------------------------------")
print()

punc = parameter_uncertainty(xopt)
for i in range(len(xopt)):
	print('{0:20s}: {1:13.5e} +- {2:13.5e}'.format(param_names[i],xopt[i],punc[i]))

interactive_plot(xopt)

















