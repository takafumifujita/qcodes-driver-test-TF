# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:49:39 2016

@author: TUD205099
"""
#%% load modules
import triple_dot
import numpy as np
import qtt
from qtt.scans import makeDataset_sweep, makeDataset_sweep_2D
import qcodes
from qcodes import load_data
#from qcodes.plots.pyqtgraph import QtPlot
from qcodes.plots.qcmatplotlib import MatPlot
from qtt.live_plotting import livePlot, fpgaCallback_2d

#%% set directory for data saving
datadir = r'K:\ns\qt\spin-qubits\data\b057_data\2017 3dot Automation\data'
qcodes.DataSet.default_io = qcodes.DiskIO(datadir)
qcodes.DataSet.default_formatter = qcodes.data.gnuplot_format.GNUPlotFormat()

#%% initialize station
remote = False

if __name__=='__main__':
    server_name = None
    station = triple_dot.initialize(server_name=server_name)
    awg = station.awg
    fpga = station.fpga
    gates = station.gates
    RF = station.RF
    keithley1 = station.keithley1
    keithley2 = station.keithley2
    keithley3 = station.keithley3
#    siggen = station.siggen
#    helium = station.helium

#%% initialize sensing dot
import qtt.structures

if __name__=='__main__' and 1:
    ggv = ['SDL', 'SDP', 'SDR']
    sdvalv = [gates.get(ggv[0]), gates.get(ggv[1]), gates.get(ggv[2])]
    sd = qtt.structures.sensingdot_t(ggv, sdvalv, station, index=1, fpga_ch=1)
    
    if 0:
        sdval, dataset = sd.autoTune(scanrange=100)
        gates.set(sd.gg[1],sdval)    
        name = '_'.join(sai.array_id for sai in dataset.keithley1_amplitude.set_arrays)
        dataset.location = dataset.location_provider(dataset.default_io, record={'name': name})
        dataset.write()

#%% defining virtual gates
from collections import OrderedDict

#L = {'P1': 1, 'P2': .454, 'P3': .19, 'D1': 1.535, 'D2': .460, 'LS': 1.020, 'RS': .167}
#L_sweep = {'P1': 1, 'P2': .454, 'P3': .19}
#M = {'P1': .537, 'P2': 1, 'P3': 0.386, 'D1': 1.206, 'D2': .855, 'LS': .318, 'RS': .251}
#M_sweep = {'P1': .537, 'P2': 1, 'P3': 0.386}
#R = {'P1': .21, 'P2': .654, 'P3': 1, 'D1': .526, 'D2': 1.493, 'LS': .164, 'RS': 1.136}
#R_sweep = {'P1': .21, 'P2': .654, 'P3': 1}

## result of capacitance measurement
# new tuning at 000-111
#L = OrderedDict([('P1', 1), ('P2', .481), ('P3', .194), ('D1', 1.224), ('D2', .380), ('LS', 0.843), ('RS', .110)])
#M = OrderedDict([('P1', .560), ('P2', 1), ('P3', .414), ('D1', 1.295), ('D2', .942), ('LS', .330), ('RS', .269)])
#R = OrderedDict([('P1', .208), ('P2', .503), ('P3', 1.25), ('D1', .394), ('D2', 1.189), ('LS', .073), ('RS', 0.914)])
L = OrderedDict([('P1', 1), ('P2', .481), ('P3', .194), ('D1', 1.224), ('D2', .380), ('LS', 0.843), ('RS', .110)])
M = OrderedDict([('P1', .540), ('P2', 1), ('P3', .414), ('D1', 1.295), ('D2', .942), ('LS', .330), ('RS', .269)])
R = OrderedDict([('P1', .208), ('P2', .503), ('P3', 1.25), ('D1', .394), ('D2', 1.189), ('LS', .073), ('RS', 0.914)])
t1 = OrderedDict([('P1', 0), ('P2', 0), ('P3', 0), ('D1', 1), ('D2', 0), ('LS', 0), ('RS', 0)])
t2 = OrderedDict([('P1', 0), ('P2', 0), ('P3', 0), ('D1', 0), ('D2', 1), ('LS', 0), ('RS', 0)])
t_L = OrderedDict([('P1', 0), ('P2', 0), ('P3', 0), ('D1', 0), ('D2', 0), ('LS', 1), ('RS', 0)])
t_R = OrderedDict([('P1', 0), ('P2', 0), ('P3', 0), ('D1', 0), ('D2', 0), ('LS', 0), ('RS', 1)])

# extraction for sweepable gates
L_sweep = OrderedDict([(list(L.keys())[0], list(L.values())[0]), (list(L.keys())[1], list(L.values())[1]), (list(L.keys())[2], list(L.values())[2])])
M_sweep = OrderedDict([(list(M.keys())[0], list(M.values())[0]), (list(M.keys())[1], list(M.values())[1]), (list(M.keys())[2], list(M.values())[2])])
R_sweep = OrderedDict([(list(R.keys())[0], list(R.values())[0]), (list(R.keys())[1], list(R.values())[1]), (list(R.keys())[2], list(R.values())[2])])

# make the cross-capacitance matrices and invert
cc_sweep = np.array([list(L_sweep.values()), list(M_sweep.values()), list(R_sweep.values())])
cc_sweep_inv = np.linalg.inv(cc_sweep)

cc = np.array([list(L.values()), list(M.values()), list(R.values()), list(t1.values()), list(t2.values()), list(t_L.values()), list(t_R.values())])
cc_inv = np.linalg.inv(cc)

# make the inverted chemical potentials
mu_L_inv = dict()
mu_M_inv = dict()
mu_R_inv = dict()
for i in range(3):
    mu_L_inv[list(L_sweep.keys())[i]] = np.dot(cc_sweep_inv, np.array([1,0,0]))[i]
    mu_M_inv[list(M_sweep.keys())[i]] = np.dot(cc_sweep_inv, np.array([0,1,0]))[i]
    mu_R_inv[list(R_sweep.keys())[i]] = np.dot(cc_sweep_inv, np.array([0,0,1]))[i]

# virtual P1+P2 gate
Dot_LR = {'P1': mu_L_inv['P1'] + mu_M_inv['P1'],
          'P2': mu_L_inv['P2'] + mu_M_inv['P2'],
          'P3': mu_L_inv['P3'] + mu_M_inv['P3']}


# dot gate dictionary for epsilon, delta, and mu
Dot_epsilon = dict()
Dot_delta   = dict()
Dot_mu      = dict()
for i in list(mu_L_inv):
    Dot_epsilon[i] = (  mu_L_inv.get(i) - mu_R_inv.get(i) ) / mu_L_inv.get('P1')
    Dot_delta[i]   = (- mu_L_inv.get(i) + mu_M_inv.get(i) - mu_R_inv.get(i) ) / mu_M_inv.get('P2')
    Dot_mu[i]      = (  mu_L_inv.get(i) + mu_M_inv.get(i) + mu_R_inv.get(i) ) / mu_L_inv.get('P1')

# dot gate dictionary for the other parameters
Dot_t1  = {list(t1.keys())[3]: cc_inv[3,3]}
Dot_t2  = {list(t2.keys())[4]: cc_inv[4,4]}
Dot_t_L = {list(t_L.keys())[5]: cc_inv[5,5]}
Dot_t_R = {list(t_R.keys())[6]: cc_inv[6,6]}
for i in range(3):
    Dot_t1[list(L_sweep.keys())[i]] = cc_inv[:,3][i]
    Dot_t2[list(L_sweep.keys())[i]] = cc_inv[:,4][i]
    Dot_t_L[list(L_sweep.keys())[i]] = cc_inv[:,5][i]
    Dot_t_R[list(L_sweep.keys())[i]] = cc_inv[:,6][i]


#%% Do 1D scan for a polarization line
fig = 1001
gg = ['P1','P3']
gatevals = gates.allvalues()
activegates = ['P1','P2','P3','D1','D2','LS','RS','SDP','SDR','SDL','T']
gate =  gg[0]
sweeprange = 5
period = 1e-3
gates.set(gg[0], 17) 
gates.set(gg[1], -211)
Naverage = 1000
waveform, sweep_info = station.awg.sweep_gate(gate, sweeprange, period)

ReadDevice = ['FPGA_ch%d' % fpga_ch]
_,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)

station.awg.stop()

dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
data = station.awg.sweep_process(dataread, waveform, Naverage)
dataset, plot = makeDataset_sweep(data, gate, sweeprange, fig=fig, gates=gates)
titletxt = plot.title.get_text()    
plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)

gates.resetgates(activegates,gatevals)
#%% Do 1D scan for a charge addition line
fig = 1002
gg = ['P2']
gate =  gg[0]
sweeprange = 5
period = 1e-3
Naverage = 1000
waittime = Naverage * period
waveform, sweep_info = station.awg.sweep_gate(gate, sweeprange, period)

ReadDevice = ['FPGA_ch%d' % fpga_ch]
_,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)

station.awg.stop()

dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
data = station.awg.sweep_process(dataread, waveform, Naverage)
dataset, plot = makeDataset_sweep(data, gate, sweeprange, fig=fig, gates=gates)
titletxt = plot.title.get_text()    
plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)
#%% Retune single gate

gates.set_P2(gates.get_P2()-0.1)

#%% Record time trace with FPGA
fig = 1005
import matplotlib.pyplot as plt
plt.close(fig)
gg = ['P2']
gate =  gg[0]
sweeprange = 0
period = 8e-3
Naverage = 1
waittime = 0
waveform, sweep_info = station.awg.sweep_gate(gate, sweeprange, period)

ReadDevice = ['FPGA_ch%d' % fpga_ch]
_,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)

station.awg.stop()

dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
data = station.awg.sweep_process(dataread, waveform, Naverage)
plot = MatPlot(data, interval=0)


#%%
from qtt.algorithms.tunneling import fit_pol_all, polmod_all_2slopes

par_fit = fit_pol_all(dataset.P1.ndarray, dataset.measured.ndarray)

#TODO: add fit to dataplot
MatPlot(dataset.P1.ndarray, polmod_all_2slopes(dataset.P1.ndarray, par_fit), interval=0)

# convert t1 from mV to GHz
t1_hz = par_fit[0]*80*(cc[0,0]-cc[1,0])/4.2

#%% focus at transition
gates.P1.set(0)
gates.P3.set(-190)

# 
gates.P1.set(17) 
gates.P3.set(-211)

#%% single fast 2D scan with virtual plungers
# TODO(TF): think how to generalize the step axis for other instruments, also for 1D
import qtt.scans
from imp import reload
reload(qtt.scans)
from qtt.scans import scan2Dfast

#scangates = ['P2','P3']
scangates = ['P1','P2','P3']
gatevals = gates.allvalues()
gg = [gates.get(scangates[0]), gates.get(scangates[1])]
activegates = ['P1','P2','P3','D1','D2','LS','RS','SDP','SDR','SDL','T']
#activegates = scangates

if __name__ == '__main__' and 1:
    stepgate = getattr(gates,scangates[0])
    stepgateval = stepgate.get()
    plot = MatPlot(interval=0)
    delay = 0.1
    scanjob = dict({'sweepdata': dict({'gate': scangates[1], 'start': gg[1] - 40, 'end': gg[1] + 40, 'step': 2.}), 'delay': delay})
    scanjob['stepdata'] = dict({'gate': scangates[0], 'start': gg[0] + 80, 'end': gg[0] - 80, 'step': -2})
    scanjob['sd'] = sd
    scanjob['sweepdata']['period'] = .5e-3
#    scanjob['gates_horz'] = {'P2':1, 'P3': 0}
#    scanjob['gates_vert'] = {'P2':0, 'P3': 1}
    scanjob['gates_horz'] = mu_R_inv
    scanjob['gates_horz'] = {'P1': 0.0066602222076373452*2, 'P2': -0.40117006534666505*2, 'P3': 0.96032257332014703*2} # mu_R_inv * const
#    scanjob['gates_vert'] = mu_M_inv
#    scanjob['gates_horz'] = Dot_epsilon
#    scanjob['gates_vert'] = Dot_mu
#    scanjob['gates_horz'] = {'P1':1, 'P2': -0.8, 'P3': -0.2}
#    scanjob['gates_vert'] = {'P1':0, 'P2': 1}
#    scanjob['gates_horz'] = {'P2':1, 'P3': 0}
#    scanjob['gates_vert'] = {'P2':0, 'P1': 1}
#    scanjob['gates_horz'] = {'P1':1, 'P3': -1.25}
#    scanjob['gates_vert'] = {'P1':1, 'P2':0.6, 'P3':1,'SDP': -.5}
#    scanjob['gates_vert'] = {'P1':1, 'P2':0.5, 'P3':1}
#    scanjob['gates_horz'] = {'P1': 1, 'P3': -1}
#    scanjob['gates_vert'] = {'P1': .5, 'P2': .5, 'P3': .5}
    scanjob['gates_vert'] = Dot_t_R
    scanjob['fpga_samp_freq'] = fpga.get_sampling_frequency()
    diff_dir = 'xy'

    RF.on()
    
    alldata = scan2Dfast(station, scanjob, liveplotwindow=plot, diff_dir=diff_dir, wait_time=None, background=False)
    plot.fig.axes[0].autoscale(tight=True)
    plot.fig.axes[1].autoscale(tight=True)
    gates.resetgates(activegates,gatevals)

    RF.off()

    if 0:
        diff_dir = 'y'
        imx = qtt.diffImageSmooth(alldata.measured.ndarray, dy=diff_dir)
        data_arr = qcodes.DataArray(name='diff', label='diff', array_id='diff', set_arrays=alldata.measured.set_arrays, preset_data=imx)
        alldata.add_array(data_arr)
    
    plot_2 = plot # reserving this plot for later analysis

#%% ADDED(TF): multiple fast-2D scans with single virtual plunger sweep, mainly for capacitance measurements
import qtt.scans
from imp import reload
reload(qtt.scans)
from qtt.scans import scan2Dfast
from qtt.tools import mouseClick

## set gate voltages to the center of the relevant charging line before running

#dotgate = 'P3' # gate corresponding to the relevant dot (fast gate)
#stepgates = ['P1','P2','D1','D2','LS','RS','SDP','SDL','SDR','T'] # all the gates for measuring cross-capacitances
#cc_init = [0.198, 0.507, 0.395, 1.184, 0.084, 0.941, 0.074, 0.064, 0.066, 3.632] # initial guess of cross capacitances

#dotgate = 'P2' # gate corresponding to the relevant dot (fast gate)
#stepgates = ['P1','P3','D1','D2','LS','RS','SDP','SDL','SDR','T'] # all the gates for measuring cross-capacitances
#cc_init = [0.537, 0.386, 1.206, 0.855, 0.318, 0.251, 0.084, 0.088, 0.066, 4.022] # initial guess of cross capacitances

dotgate = 'P1' # gate corresponding to the relevant dot (fast gate)
stepgates = ['P2','P3','D1','D2','LS','RS','SDP','SDL','SDR','T'] # all the gates for measuring cross-capacitances
cc_init = [0.454, 0.190, 1.535, 0.460, 1.020, 0.167, 0.079, 0.065, 0.067, 3.186] # initial guess of cross capacitances

try:
    len(stepgates) == len(cc_init)
except ValueError:
    print("Oops!  That was no valid VECTOR.  Try again...")

RF.on()
for step_num, stepgate_name in enumerate(stepgates):
    scangates = [stepgate_name,dotgate]
    gatevals = gates.allvalues()
    gg = [gates.get(scangates[0]), gates.get(scangates[1])]
    activegates = scangates

    if __name__ == '__main__' and 1:
        stepgate = getattr(gates,scangates[0])
        stepgateval = stepgate.get()
        plot = MatPlot(interval=0)
        delay = 0.1
        step_width = 12 / 2 
        scanjob = dict({'sweepdata': dict({'gate': scangates[1], 'start': gg[1] - step_width*cc_init[step_num], 'end': gg[1] + step_width*cc_init[step_num], 'step': 1.}), 'delay': delay})
        scanjob['stepdata'] = dict({'gate': scangates[0], 'start': gg[0] + step_width, 'end': gg[0] - step_width, 'step': -0.1})
        scanjob['sd'] = sd
        scanjob['sweepdata']['period'] = .5e-3
        if dotgate == 'P2':
            scanjob['gates_horz'] = {dotgate:1, 'P1': 0} # second one is a dummy gate
        else:
            scanjob['gates_horz'] = {dotgate:1, 'P2': 0} # second one is a dummy gate
        scanjob['gates_vert'] = {dotgate:0, stepgate_name: 1}
        scanjob['fpga_samp_freq'] = fpga.get_sampling_frequency()
        diff_dir = 'xy'
    
        alldata = scan2Dfast(station, scanjob, liveplotwindow=plot, wait_time=None, background=False)
        plot.fig.axes[0].autoscale(tight=True)
        plot.fig.axes[1].autoscale(tight=True)
        gates.resetgates(activegates,gatevals)
    
        exec('plot_' + str(step_num) + '=plot') # reserving plot for later analysis
        clicks = mouseClick(plot) # run mauseClick() class from the develp_automated_crosscapacitance.py
RF.off()

#%%
for g in scanjob['gates_vert']:
    gates.get(g)
    
#%% do a single qcodes loop of the sensing dot plunger
if __name__=='__main__' and 1:
    SDP_val = gates.SDP.get()
    RF.on()
    station.set_measurement(station.keithley1.amplitude)
    loop_1d = qcodes.Loop(gates.SDP[-250:-300:1],delay=0.1)
    dataset = loop_1d.run(background=False, data_manager=False)
    qcodes.plots.qcmatplotlib.MatPlot(dataset.default_parameter_array(), interval=0)
    gates.SDP.set(SDP_val)

#%% 2d scan
if __name__=='__main__' and 1:
#    loop_2d = qcodes.Loop(gates.SDL[-300:0:5],delay=0.1).loop(gates.SDR[-10:0:5])
    loop_2d = qcodes.Loop(gates.SDL[-300:0:5],delay=0.1).each(qcodes.Loop(gates.SDR[-300:0:5],delay=.1).each(station.keithley2.amplitude))
    dataset_2d = loop_2d.run(background=False, data_manager=False)
    qcodes.plots.qcmatplotlib.MatPlot(dataset_2d.default_parameter_array(), interval=0)
    
#%% single veryfast 2D scan
if __name__=='__main__':
    fig = 111
    sweepgates = ['P1','P3']
    sweepranges = [80, 80]
    resolution = [90,90]
    Naverage = 1000
    fpga_ch = 1
    diff_dir = 'xy'
#    diff_dir = None
    
    waveform, sweep_info = station.awg.sweep_2D(station.fpga.get_sampling_frequency(), sweepgates, sweepranges, resolution)
    waittime = resolution[0]*resolution[1]*Naverage/fpga.get_sampling_frequency()
    
    ReadDevice = ['FPGA_ch%d' % fpga_ch]
    _,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)
    
    station.awg.stop()
    
    dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
    data = station.awg.sweep_2D_process(dataread, waveform, diff_dir=diff_dir)
    dataset, plot = makeDataset_sweep_2D(data, gates, sweepgates, sweepranges, fig=111)
    titletxt = plot.title.get_text()    
    plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)
    
    name = '_'.join(sai.array_id for sai in dataset.measured.set_arrays)
    dataset.location = dataset.location_provider(dataset.default_io, record={'name': name})
    dataset.write()

#%%
station.awg.sweep_run(sweep_info)
station.awg.stop()

dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
data = station.awg.sweep_2D_process(dataread, waveform, diff_dir=diff_dir)
dataset, plot = makeDataset_sweep_2D(data, gates, sweepgates, sweepranges, fig=111)
titletxt = plot.title.get_text()    
plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)

name = '_'.join(sai.array_id for sai in dataset.measured.set_arrays)
dataset.location = dataset.location_provider(dataset.default_io, record={'name': name})
dataset.write()

#%% single veryfast 2D scan of virtual gates
if __name__=='__main__':
    fig = 203
#    gates_horz = {'P1': 1, 'P3': -.0}
#    gates_vert = {'P3': 1, 'P1': -.0}
#    gates_horz = {'P1': 1, 'P3': -1}
#    gates_vert = {'P1': -.5, 'P2': 1, 'P3': -.5}
#    gates_horz = {'P1': 1, 'P3': -1}
#    gates_vert = {'P1': .5, 'P2': .5, 'P3': .5}
#    gates_horz = mu_M_inv
#    gates_vert = mu_R_inv
    gates_horz = Dot_epsilon
    gates_vert = Dot_delta
    gates_vert = Dot_mu
    sweepranges = [90, 150]
    sweepgates = ['P1', 'P2']
    resolution = [90,90]
    Naverage = 1000
    fpga_ch = 1
    diff_dir = 'xy'
#    diff_dir = None    
    
    waveform, sweep_info = station.awg.sweep_2D_virt(station.fpga.get_sampling_frequency(), gates_horz, gates_vert, sweepranges, resolution)
    waittime = resolution[0]*resolution[1]*Naverage/fpga.get_sampling_frequency()
    
    ReadDevice = ['FPGA_ch%d' % fpga_ch]
    _,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)
    
    station.awg.stop()
    
    dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
    data = station.awg.sweep_2D_process(dataread, waveform, diff_dir=None)
    dataset, plot = makeDataset_sweep_2D(data, gates, sweepgates, sweepranges, fig=fig)
    titletxt = plot.title.get_text()    
    plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)
    plot.fig.axes[0].autoscale(tight=True)
    plot.fig.axes[1].autoscale(tight=True)
    
    fig = 204
    dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
    data = station.awg.sweep_2D_process(dataread, waveform, diff_dir=diff_dir)
    dataset, plot = makeDataset_sweep_2D(data, gates, sweepgates, sweepranges, fig=fig)
    titletxt = plot.title.get_text()    
    plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)
    plot.fig.axes[0].autoscale(tight=True)
    plot.fig.axes[1].autoscale(tight=True)
    
    name = '_'.join(sai.array_id for sai in dataset.measured.set_arrays)
    dataset.location = dataset.location_provider(dataset.default_io, record={'name': name})
    dataset.write()

#%% videomode tuning
if __name__=='__main__':
    sweepgates = ['P1','P3']
    sweepranges = [100,100]
    resolution = [90,90]
    Naverage = 25
    fpga_ch = 1
    diff_dir = 'xy'
    
    waveform, sweep_info = station.awg.sweep_2D(station.fpga.get_sampling_frequency(), sweepgates, sweepranges, resolution)

    lp = livePlot(gates, sweepgates, sweepranges)
    lp.datafunction = fpgaCallback_2d(station, waveform, Naverage, fpga_ch, resolution, diff_dir)
    lp.startreadout(rate=10)

#%% videomode tuning with virtual gates
if __name__=='__main__':
    sweepgates = ['P1','P3']
#    gates_horz = {'P1':1, 'P3': -1.25}
#    gates_vert = {'P1':1, 'P2':.5, 'P3':0.8}
##    gates_horz = {'P1': 1, 'P3': -.2}
##    gates_vert = {'P3': 1, 'P1': -.2}
#    gates_horz = {'P1': 1, 'P3': -1.25}
#    gates_vert = {'P1': -.5, 'P2': 1, 'P3': -.4}
#    gates_horz = {'P1': 1, 'P2': 0}
#    gates_vert = {'P1': 0, 'P3': 1}
##    gates_horz = mu_R_inv
    gates_horz = mu_R_inv
    gates_vert = Dot_LR
#    gates_horz = Dot_epsilon
#    gates_vert = Dot_mu
    sweepranges = [90, 140]
    sweepranges = [140, 90]
    resolution = [90,90]
    Naverage = 25
    fpga_ch = 1
    diff_dir = 'xy'
    RF.on()

    waveform, sweep_info = station.awg.sweep_2D_virt(station.fpga.get_sampling_frequency(), gates_horz, gates_vert, sweepranges, resolution)

    lp = livePlot(gates, sweepgates, sweepranges)
    lp.datafunction = fpgaCallback_2d(station, waveform, Naverage, fpga_ch, resolution, diff_dir)
    lp.startreadout(rate=10)
    
#%% stop video mode
lp.stopreadout()
awg.stop()
RF.off()

#%% single 2D scan (2D scan in scans does not work)
# crashes still! (does it matter which direction we step?)
scangates = ['SDL','SDR']
gg = [0, 0]
#gg = [136, 100]
#gg = [gates.get(sweepgates[0]), gates.get(sweepgates[1])]

if __name__ == '__main__' and 1:
    gate_horz = getattr(gates, scangates[0])
    gate_vert = getattr(gates, scangates[1])
    delay = .01
    scanjob = dict({'sweepdata': dict({'gate': scangates[1], 'start': gg[1] - 500, 'end': gg[1] , 'step': 5.}), 'delay': delay})
    scanjob['stepdata'] = dict({'gate': scangates[0], 'start': gg[0] - 500, 'end': gg[0], 'step': 5.})    
    
    # move/combine this code to scan2D in qtt.scans
    station.set_measurement(station.keithley2.amplitude)
    loop2D = qcodes.Loop(gate_horz[scanjob['sweepdata']['start']:scanjob['sweepdata']['end']:scanjob['sweepdata']['step']], delay=delay)
    loop2D_full = loop2D.loop(gate_vert[scanjob['stepdata']['start']:scanjob['stepdata']['end']:scanjob['stepdata']['step']], delay=delay)
    alldata = loop2D_full.run(background=False, data_manager=False)

    gates.set(scangates[0], gg[0])
    gates.set(scangates[1], gg[1])

#%% scan all pairs of plungers for dots below
plungers = ['P1','P2','P3']
fpga.set_sampling_frequency(200000)

pairs = []
for i in plungers:
    for j in plungers:
        if i != j:
            pairs += [[i,j]]

if __name__=='__main__':
    plungers = ['P1','P2','P3']
    fignum = 500
    for i in range(0,len(plungers)*2):
        fig = fignum + i
        sweepgates = pairs[i]
        sweepranges = [80, 80]
        resolution = [90,90]
        Naverage = 1000
        fpga_ch = 1
        diff_dir = 'xy'
        
        waveform, sweep_info = station.awg.sweep_2D(station.fpga.get_sampling_frequency(), sweepgates, sweepranges, resolution)
        waittime = resolution[0]*resolution[1]*Naverage/fpga.get_sampling_frequency()
        
        ReadDevice = ['FPGA_ch%d' % fpga_ch]
        _,DataRead_ch1,DataRead_ch2 = station.fpga.readFPGA(Naverage=Naverage, ReadDevice=ReadDevice, waittime=waittime)
        
        station.awg.stop()
        
        dataread = [DataRead_ch1,DataRead_ch2][fpga_ch-1]
        data = station.awg.sweep_2D_process(dataread, waveform, diff_dir=diff_dir)
        dataset, plot = makeDataset_sweep_2D(data, gates, sweepgates, sweepranges, fig=fig)
        titletxt = plot.title.get_text()    
        plot.title.set_text(titletxt + ', diff_dir: %s' % diff_dir)
        
        name = '_'.join(sai.array_id for sai in dataset.measured.set_arrays)
        dataset.location = dataset.location_provider(dataset.default_io, record={'name': name})
        dataset.write()
    
    qtt.pmatlab.tilefigs([500,501,502,503,504,505],[3,2])

#%% Fit lines and calculate slopes
from qtt.deprecated.linetools import costFunctionLine
from scipy.optimize import minimize
import matplotlib.pyplot as plt

pp = ['P3', 'P1']
verbose = 0
cb = None
param0 = [-5,5,.5*np.pi] # x,y,theta, 
px = [dataset.measured.ndarray.shape[0]//2,dataset.measured.ndarray.shape[1]//2]
imx = dataset.measured.ndarray
istep = .5
cgate = pp[0]
igate = pp[1]

costfun = lambda x : costFunctionLine(x, -imx, verbose=0, istep=istep, px=px, dthr=1, dwidth=2)
res = minimize(costfun, param0, method='powell', options={'maxiter': 3000, 'maxfev': 101400, 'xtol': 1e-8, 'disp': verbose>=2}, callback=cb)
c = costFunctionLine(res.x, imx, istep, verbose=1, fig=fig, px=px); plt.figure(fig); plt.xlabel(cgate); plt.ylabel(igate); plt.close(fig+1)

#%% parameterviewer
from qtt.parameterviewer import createParameterWidgetRemote
from qtt.parameterviewer import createParameterWidget

if __name__=='__main__' and not remote:
    p = createParameterWidget([gates,])

if __name__=='__main__' and remote:
    p=createParameterWidgetRemote([gates,])


#%% load data and plot results
if __name__=='__main__' and 1:
    olddatadir = r'K:\ns\qt\spin-qubits\data\b057_data\2016 3dot experiment\data\2016-11-11\18-20-44_P2_P3'
    dataset_old = load_data(location=olddatadir)
#    qcodes.plots.pyqtgraph.QtPlot(dataset_old.measured, interval=0)
    plotje = qcodes.plots.qcmatplotlib.MatPlot(dataset_old.measured, interval=0)

#%% delta
x = 2
gates.P1.set(gates.P1.get()+Dot_delta['P1']*x)
gates.P2.set(gates.P2.get()+Dot_delta['P2']*x)
gates.P3.set(gates.P3.get()+Dot_delta['P3']*x)

#%% epsilon
x = 2
gates.P1.set(gates.P1.get()+Dot_epsilon['P1']*x)
gates.P2.set(gates.P2.get()+Dot_epsilon['P2']*x)
gates.P3.set(gates.P3.get()+Dot_epsilon['P3']*x)

#%% mu
x = -10
gates.P1.set(gates.P1.get()+Dot_mu['P1']*x)
gates.P2.set(gates.P2.get()+Dot_mu['P2']*x)
gates.P3.set(gates.P3.get()+Dot_mu['P3']*x)
gates.SDP.set(gates.SDP.get() - 0.15*x)


#%% tuning t1
x = 40
gates.D1.set(gates.D1.get()+Dot_t1['D1']*x)
gates.P1.set(gates.P1.get()+Dot_t1['P1']*x)
gates.P2.set(gates.P2.get()+Dot_t1['P2']*x)
gates.P3.set(gates.P3.get()+Dot_t1['P3']*x)
gates.allvalues()

#%% tuning t2
x = 10
gates.D2.set(gates.D2.get()+Dot_t2['D2']*x)
gates.P1.set(gates.P1.get()+Dot_t2['P1']*x)
gates.P2.set(gates.P2.get()+Dot_t2['P2']*x)
gates.P3.set(gates.P3.get()+Dot_t2['P3']*x)
gates.allvalues()

#%% tuning t_L
x = -20
gates.LS.set(gates.LS.get()+Dot_t_L['LS']*x)
gates.P1.set(gates.P1.get()+Dot_t_L['P1']*x)
gates.P2.set(gates.P2.get()+Dot_t_L['P2']*x)
gates.P3.set(gates.P3.get()+Dot_t_L['P3']*x)
gates.allvalues()

#%% tuning t_R
x = 80
gates.RS.set(gates.RS.get()+Dot_t_R['RS']*x)
gates.P1.set(gates.P1.get()+Dot_t_R['P1']*x)
gates.P2.set(gates.P2.get()+Dot_t_R['P2']*x)
gates.P3.set(gates.P3.get()+Dot_t_R['P3']*x)
gates.allvalues()

#%% L
x = -5
gates.P1.set(gates.P1.get()+mu_L_inv['P1']*x)
gates.P2.set(gates.P2.get()+mu_L_inv['P2']*x)
gates.P3.set(gates.P3.get()+mu_L_inv['P3']*x)

#%% M
x = -10
gates.P1.set(gates.P1.get()+mu_M_inv['P1']*x)
gates.P2.set(gates.P2.get()+mu_M_inv['P2']*x)
gates.P3.set(gates.P3.get()+mu_M_inv['P3']*x)

#%% R
x = 30
gates.P1.set(gates.P1.get()+mu_R_inv['P1']*x)
gates.P2.set(gates.P2.get()+mu_R_inv['P2']*x)
gates.P3.set(gates.P3.get()+mu_R_inv['P3']*x)

#%% sensing dot
x = -3
gates.SDP.set(gates.SDP.get() + x)

#%% LS,RS?
x = 5
gates.P1.set(gates.P1.get()-x)
gates.P2.set(gates.P2.get()+.6*x)
gates.P3.set(gates.P3.get()-x)
gates.LS.set(gates.LS.get()+x)
gates.RS.set(gates.RS.get()+x)

#%%
x = 5
gates.P1.set(gates.P1.get()+x*1)
gates.P2.set(gates.P2.get()+x*.5)
gates.P3.set(gates.P3.get()+x*.8)

#%%
x = 10
gates.T.set(gates.T.get()+x)
gates.P1.set(gates.P1.get()-x)
gates.P2.set(gates.P2.get()-x)
gates.P3.set(gates.P3.get()-x)
gates.D1.set(gates.D1.get()-x/2)
gates.D2.set(gates.D2.get()-x/2)
gates.LS.set(gates.LS.get()-x/2)
gates.RS.set(gates.RS.get()-x/2)
gates.SDR.set(gates.SDR.get()-x)
gates.SDP.set(gates.SDP.get()-x)
gates.SDL.set(gates.SDL.get()-x)

#%%
x = -10
gates.D2.set(gates.D2.get()+x)
gates.P2.set(gates.P2.get()-x/2)
gates.P3.set(gates.P3.get()-x/2)

#%%
# pid 5604: WARNING helpers.py:204 - negative delay -0.054258 sec

#%% FPGA marker error?
#Sending the waveform sweep_P1
#Sending the waveform sweep_P3
#scan2Dfast: 0/50: setting P3 to -189.059
#Traceback (most recent call last):
#
#  File "<ipython-input-6-b0f69151b5a6>", line 40, in <module>
#    alldata = scan2Dfast(station, scanjob, liveplotwindow=plot, diff_dir=diff_dir, wait_time=None, background=False)
#
#  File "D:\Users\diepencjv\qtt\qtt\scans.py", line 492, in scan2Dfast
#    alldata.measured.ndarray[ix] = readfunc(waveform, Naverage)
#
#ValueError: could not broadcast input array from shape (462) into shape (16)

#%% voltages to reset to
basevalues = {'D1': 0.030518043793335892,
 'D2': 0.030518043793335892,
 'LS': 0.030518043793335892,
 'LS_fine': 0.030518043793335892,
 'P1': 0.030518043793335892,
 'P1_fine': 0.21362630655380599,
 'P2': 0.030518043793335892,
 'P2_fine': 1.0070954451819034,
 'P3': 0.030518043793335892,
 'P3_fine': 4.242008087281647,
 'QPC': 0.030518043793335892,
 'RS': 0.030518043793335892,
 'RS_fine': 0.030518043793335892,
 'SDL': -326.02426184481578,
 'SDP': 0.030518043793335892,
 'SDP_fine': 0.030518043793335892,
 'SDR': -389.99008163576718,
 'T': -95.796139467460307,
 'bias_1': 0.030518043793335892,
 'bias_2': -499.97711146715505,
 'bias_3': 0.030518043793335892,
 'bias_4': 0.030518043793335892}



#%% reset gates to basevalues
activegates = basevalues.keys()
gates.resetgates(activegates, basevalues)