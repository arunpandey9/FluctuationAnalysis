# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:29:41 2020

@author: apan
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:04:51 2019

@author: luru
"""

from W7X_TDU_LP_access import archivedb_access as adba
#from bdbs_MDSplus_W7X import faster_time_context as mdspdba

import numpy as np
import copy
import time
import os
import matplotlib.pyplot as plt

from char_compare import utilities as ut
from char_compare import probe_geometry as pg
from char_compare import plotting_utilities as put
from char_compare import Probe
# from char_compare import test_data as td


class Manager:
    def __init__(self, program_no, time_range, **kwargs):
        self._program_no = program_no
        self.time_range = time_range
        if not hasattr(self, '_probes'):
            self._probes = {}

    def __new__(cls, program_no, time_range, **kwargs):
        from_cache = kwargs.pop('from_cache', False)
        path = kwargs.pop('path', None)
        info = (program_no, time_range, kwargs.pop('content', 'all'))
        if from_cache:
            try:
                return cls.load(path=path, info=info)
            except FileNotFoundError:
                print('File not found, loading from DB')
        return super(Manager, cls).__new__(cls)

    def __getnewargs__(self):
        return self.program_number, self.time_range

    @property
    def program_number(self):
        return self._program_no

    @property
    def time_range(self):
        return self._time_range

    @time_range.setter
    def time_range(self, time_range):
        if len(time_range) != 2:
            raise NotImplementedError()
        time_range = [float(t) for t in time_range]
        self._time_range = time_range

    def getitem(self, comp_number_tuple, **kwargs):
        """Nachtraegliches zuladen geht nicht.

        kwargs:
            verbose (int)
            force_supply_voltages (boolean)
            load_raw (boolean)
            load_processed (boolean | list(str))
        """
        verbose = kwargs.get('verbose', 0)
        force_supply_voltages = kwargs.pop('force_supply_voltages', False)
        load_raw = kwargs.pop('load_raw', True)
        load_processed = kwargs.pop('load_processed',
                                    ['Te', 'ne', 'beta', 'Vf'])
        if (load_processed is not None and isinstance(load_processed, str)):
            load_processed = [load_processed]
        component, number = comp_number_tuple
        if (component, number) not in self.probes:
            new_probe = Probe(component, number, None, None, None,
                                self.program_number)
            if not (load_processed or load_raw):
                ut.vprint('Loading nothing', verbose, activation=0)
            if load_raw:
                ut.vprint('Loading raw data', verbose)
                if component == 'test':
                    # get testdata
                    raise NotImplementedError()
                    # return_dictionary = td.generate_test_data(number)
                elif component in ('MLP', 'BLP'):
                    raise NotImplementedError
                    # get MLP or BLP signal
                    '''
                    return_dictionary = {}
                    for q in 'IU':
                        ADC = int(input('Input ADC number for {}: '.format(q)))
                        Channel = int(
                            input('Input channel number for {}: '.format(q)))
                        time, data = mdspdba.getsegdata(
                            shot=self.program_number,
                            treename='QRP',
                            nodename='HARDWARE.SIS8300KU_{}.CHANNEL{}'.format(
                                ADC, Channel),
                            time_range=self.time_range,
                            verbose=0)
                        return_dictionary['{}_{}'.format(q, number)] = data
                        return_dictionary['{}_t_{}'.format(q, number)] = time
                    '''
                else:
                    # get data for specific CLP from archive
                    return_dictionary = adba.retrieve(
                        self.program_number,
                        component,
                        probes=[number],
                        signals=['I', 'U'],
                        time_range=self.time_range, verbose=verbose)

                    if (('U_{}'.format(number) not in return_dictionary)
                            or force_supply_voltages):
                        # get power supply voltage if it wasn't directly measured
                        power_supply_name = adba.find_power_supply(
                            component, number, self.program_number)
                        power_supply_voltage = adba.retrieve(
                            self.program_number,
                            'outsideVessel',
                            probes=[power_supply_name],
                            signals=['U'],
                            time_range=self.time_range)['U_' +
                                                        power_supply_name]
                        return_dictionary['U_{}'.format(
                            number)] = power_supply_voltage

                min_length = np.min(
                    [len(v) for v in return_dictionary.values()])

                if 'I_{}'.format(number) not in return_dictionary:
                    print('No current recorded for probe {}'.format(number))
                    return_dictionary['I_{}'.format(number)] = np.zeros(1)
                    return_dictionary['I_t_{}'.format(
                        number)] = return_dictionary['U_t_{}'.format(number)]

                new_probe = Probe(
                    component, number,
                    return_dictionary['I_t_{}'.format(number)][:min_length],
                    return_dictionary['I_{}'.format(number)][:min_length],
                    return_dictionary['U_{}'.format(number)][:min_length],
                    self.program_number)

            if load_processed:
                time = []
                new_probe._fit_results['minerva', None] = {}
                new_probe._fit_results['minerva', None]['params'] = {}
                new_probe._fit_results['minerva', None]['variance'] = {}
                for q in load_processed:
                    time, new_probe._fit_results[
                        'minerva', None]['params'][q] = adba.load_processed(
                            prog=self.program_number,
                            component=component,
                            quantity=q,
                            probe=number,
                            time_range=self.time_range,
                            **kwargs)
                    time, new_probe._fit_results[
                        'minerva', None]['variance'][q] = adba.load_processed(
                            prog=self.program_number,
                            component=component,
                            quantity=q + '_std',
                            probe=number,
                            time_range=self.time_range,
                            **kwargs)
                new_probe._fit_results['minerva', None]['time'] = time
                # if new_probe.time is None: new_probe.time = time #unnecessary and unwise
            self.probes[component, number] = new_probe
        return self.probes[component, number]

    def __getitem__(self, comp_number_tuple):
        return self.getitem(comp_number_tuple)

    @property
    def probes(self):
        return self._probes

    @probes.setter
    def probes(self, probes):
        self._probes = probes

    def keys(self):
        return list(self.probes)

    def copy(self):
        return copy.copy(self)

    def keep(self, mode):
        for pi in self.probes:
            self.probes[pi].keep(mode)

    def save(self, path=None, content='all'):
        if path is None:
            directory = 'C:/Users/luru/Desktop/GitHub/Data'
            if self.program_number not in os.listdir(directory):
                os.mkdir(directory + '/' + self.program_number)
            path = directory + '/{pn}/probes_{tr}_{c}.npz'.format(
                pn=self.program_number, tr=self.time_range, c=content)
        if content != 'all':
            obj = self.copy()
            obj.keep(mode=content)
            obj.save(path)
            return  # is this ever reached?
        np.savez(path, self)

    @classmethod
    def load(cls, path=None, info=None):
        '''Specify either path or info=(program_number,time_range,mode)'''
        if path is None:
            if info is None:
                raise Exception("Don't know what to laod, please provide info")
            path = 'C:/Users/luru/Desktop/GitHub/Data/{}/probes_{}_{}.npz'.format(
                *info)
        return np.load(path)['arr_0'].tolist()

    def load_probes(self, component, probes, **kwargs):
        new_probes_flag = False
        for pi in probes:
            if (component, pi) not in self.probes:
                self.getitem((component, pi), **kwargs)
                new_probes_flag = True
        return new_probes_flag

    def do_all(self, func, *args, **kwargs):
        timed = kwargs.pop('timed', False)
        verbose = kwargs.get('verbose', 0)
        kwargs['verbose'] = verbose - 1
        mute_warnings = kwargs.pop('mute_warnings', False)
        np.warnings.filterwarnings('ignore')  # Fuck it for now
        if timed:
            start = time.time()
        for p_key, p_obj in self.probes.items():
            # args=[ut.instanciate_with_probe(a,p_obj) for a in initial_args]
            if verbose > 0:
                print('Applying {} to {}'.format(func, p_key))
            if timed:
                loop_start = time.time()
            if mute_warnings:
                with np.warnings.catch_warnings():
                    [np.warnings.filterwarnings(*m_w) for m_w in mute_warnings]
                    getattr(p_obj, func)(*args, **kwargs)
            else:
                getattr(p_obj, func)(*args, **kwargs)
            if timed:
                completion_time = time.time()
                print('Applying {} to {} took {:.2n}s, running since {:.2n}s'.
                      format(func, p_key, completion_time - loop_start,
                             completion_time - start))

    def plot_overview(self, model, slicer, component, quantity, **kwargs):
        profile_cut = kwargs.pop('profile_cut', None)
        time_cut = kwargs.pop('time_cut', None)
        method = kwargs.pop('method', 'pcolormesh')
        save = kwargs.pop('save', True)
        probes = kwargs.pop('probes', np.arange(1, 11))
        verbose = kwargs.pop('verbose', 0)
        super_q = kwargs.pop('super_q', '')

        # load everything
        for pi in probes:
            if (component, pi) not in self.probes:
                self.__getitem__((component, pi), **kwargs)
        probe_nums = self._probe_numbers(component)
        if profile_cut:
            if not hasattr(profile_cut, '__iter__'):
                profile_cut = [profile_cut]
            for i, pc in enumerate(profile_cut):
                if pc < self.time_range[0]:
                    ut.vprint('Setting profile cut to edge', verbose)
                    profile_cut[i] = self.time_range[0]
                elif pc > self.time_range[1]:
                    ut.vprint('Setting profile cut to edge', verbose)
                    profile_cut[i] = self.time_range[1]

        if time_cut:
            if not hasattr(time_cut, '__iter__'):
                time_cut = [time_cut]
            for i, tc in enumerate(time_cut):
                if not (tc in probe_nums):
                    ut.vprint('Setting time cut to highest available', verbose)
                    time_cut[i] = max(probe_nums)

        fig = plt.figure('{} Profile {} {} {}'.format(quantity, component,
                                                      self.program_number,
                                                      self.time_range))
        # if this is nasty use gridSpec instead
        # setup all cases
        if profile_cut and time_cut:
            ut.vprint('plot all', verbose)
            ax_2D = plt.subplot2grid((3, 4), (0, 0), colspan=3, rowspan=2)
            ax_cb = plt.subplot2grid((3, 4), (2, 3), colspan=1, rowspan=1)
            ax_profile = plt.subplot2grid((3, 4), (0, 3),
                                          colspan=1,
                                          rowspan=2,
                                          sharey=ax_2D)
            ax_time = plt.subplot2grid((3, 4), (2, 0),
                                       colspan=3,
                                       rowspan=1,
                                       sharex=ax_2D)

            plt.setp(ax_2D.get_xticklabels(), visible=False)
            ax_time.set_ylabel('{} [{}]'.format(quantity,
                                                ut.unit_dictionary[quantity]))
            ax_time.set_xlabel('Time [s]')
            ax_time.set_title('Timetrace for probe {}'.format(time_cut))
            ax_profile.set_title('Profile at {}'.format(profile_cut))
            ax_profile.yaxis.tick_right()
            ax_profile.yaxis.set_label_position('right')
            ax_profile.set_ylabel('Distance from pumping gap [mm]')
        elif profile_cut and not time_cut:
            ut.vprint('plot profile', verbose)
            ax_2D = plt.subplot2grid((3, 8), (0, 0), colspan=7, rowspan=2)
            ax_cb = plt.subplot2grid((3, 8), (0, 7), colspan=1, rowspan=2)
            ax_cb.tick_params(left=False, right=True)
            ax_profile = plt.subplot2grid((3, 8), (2, 0), colspan=7, rowspan=1)
            ax_profile.set_ylabel('{} [{}]'.format(
                quantity, ut.unit_dictionary[quantity]))
            ax_profile.set_title('Profile at {}'.format(profile_cut))
            ax_profile.set_xlabel('Distance from pumping gap [mm]')
        elif time_cut and not profile_cut:
            ut.vprint('plot time', verbose)
            ax_2D = plt.subplot2grid((3, 8), (0, 0), colspan=7, rowspan=2)
            ax_cb = plt.subplot2grid((3, 8), (0, 7), colspan=1, rowspan=2)
            ax_cb.tick_params(left=False, right=True)
            ax_time = plt.subplot2grid((3, 8), (2, 0),
                                       colspan=7,
                                       rowspan=1,
                                       sharex=ax_2D)
            plt.setp(ax_2D.get_xticklabels(), visible=False)
            ax_time.set_ylabel('{} [{}]'.format(quantity,
                                                ut.unit_dictionary[quantity]))
            ax_time.set_xlabel('Time [s]')
            ax_time.set_title('Timetrace for probe {}'.format(time_cut))
        else:
            ut.vprint('plot only 2D', verbose)
            ax_2D = plt.subplot2grid((3, 8), (0, 0), colspan=7, rowspan=3)
            ax_cb = plt.subplot2grid((3, 8), (0, 7), colspan=1, rowspan=3)
            ax_2D.set_xlabel('Time [s]')
        # Always
        ax_2D.set_ylabel('Distance from pumping gap [mm]')
        ax_2D.set_title('{} {} {}'.format(quantity, component,
                                          self.program_number))

        # 2D plot
        time, positions, values = self.get_quantity(model, slicer, component,
                                                    quantity, **kwargs)
        ax_t, ax_p = np.meshgrid(ut.continue_sequence(time),
                                 ut.continue_sequence(positions))
        contours = getattr(ax_2D, method)(ax_t, ax_p, values)
        # colorbar
        cb = fig.colorbar(contours, cax=ax_cb, aspect=1, drawedges=False)
        cb.set_label('{} in [{}]'.format(quantity,
                                         ut.unit_dictionary[quantity]))

        if profile_cut is not None:
            color_wheel = put.color_wheel(5)
            for pc in profile_cut:
                color = next(color_wheel)
                if time_cut:
                    ax_profile.plot(
                        values[:, np.argwhere(time >= pc)[0]].flatten(),
                        positions,
                        marker='o',
                        color=color)
                else:
                    ax_profile.plot(
                        positions,
                        values[:, np.argwhere(time >= pc)[0]].flatten(),
                        marker='o',
                        color=color)
                ax_2D.vlines(pc, ax_p[0, 0], ax_p[-1, 0], color=color)
        if time_cut is not None:
            color_wheel = put.color_wheel(5)
            for tc in time_cut:
                pos = pg.distance_to_gap(component, tc)
                color = next(color_wheel)
                ax_time.plot(
                    time,
                    values[np.argwhere(positions >= pos)[-1], :].flatten(),
                    color=color)
                ax_2D.hlines(pos, ax_t[0, 0], ax_t[0, -1], color=color)

        if save:
            # use plot_utils.savefig
            plt.gcf().set_size_inches(19, 11)
            path = 'C:/Users/luru/Desktop/GitHub/Data/{pn}/'
            # file_name='{}_profile_{}_{}.png'
            file_name = '{sq}{q}_profile_{c}_{tr}_{m}_{s}.png'
            plt.savefig((path + file_name).format(q=quantity,
                                                  pn=self.program_number,
                                                  tr=self.time_range,
                                                  m=str(model),
                                                  s=str(slicer),
                                                  c=component,
                                                  sq=super_q))

        return fig, fig.get_axes(), cb

    def get_quantity(self, model, slicer, component, quantity, **kwargs):
        # if no specific time probe is set just pick the first one from the dict
        time_probe = kwargs.pop('time_probe', 0) 
        # time_fun_kwarg = kwargs.pop('time_fun_kwarg', 'first')
        # time_fun = ut.time_fun(time_fun_kwarg)
        super_q = kwargs.pop('super_q', 'params')
        if time_probe == 0:
            time_probe = self._probe_numbers(component=component)[0]
        values = []
        for pi in self.probes:
            values.append(self[pi].fit_result(model,
                                              slicer,
                                              super_q,
                                              sub_key=quantity))
        time = self[component, time_probe].fit_result(model, slicer, 'time')
        positions = pg.distances_to_gap(component,
                                        self._probe_numbers(component))
        return np.asarray(time), np.asarray(positions), np.asarray(values)

    def _probe_numbers(self, component=None):
        names = self.probes.keys()
        if component:
            return [i for c, i in names if c == component]
        else:
            return [i for c, i in names]

    def save_quantity(self, model, slicer, component, quantity, **kwargs):
        path = kwargs.pop('path', None)
        directory = kwargs.pop('directory', None)
        super_q = kwargs.get('super_q', '')
        if directory is None:
            directory = 'C:/Users/luru/Desktop/GitHub/Data'
        if self.program_number not in os.listdir(directory):
            os.mkdir(directory + '/' + self.program_number)
        if path is None:
            path = directory + '/{pn}/{c}_{sq}{q}_{tr}_{m}_{s}.npz'.format(
                q=quantity,
                pn=self.program_number,
                tr=self.time_range,
                m=model.name,
                s=str(slicer),
                c=component,
                sq=super_q)

        np.savez(
            path,
            *self.get_quantity(model, slicer, component, quantity, **kwargs))

    def saveStuff(self, key, component, **kwargs):
        quantities = kwargs.pop('quantities',['Te','ne','Vf','beta'])
        for q in quantities:
            self.save_quantity(*key, component, q, **kwargs)    
            self.save_quantity(*key, component, q, super_q='variance', **kwargs)

    def displayStuff(self, key, component, **kwargs):
        self.plot_overview(*key, component, 'Te', **kwargs)
        self.plot_overview(*key, component, 'ne', **kwargs)
        self.plot_overview(*key, component, 'Vf', **kwargs)
        self[component, 10].plot_many_characteristics(*key)

# =============================================================================
#     def savePlots(self,key):
#         self.plot_overview(*key,ut.component_dict['LD'],'Te')
#         self.plot_overview(*key,ut.component_dict['LD'],'ne')
#         self[ut.component_dict['LD'],10].plot_many_characteristics(*key)
# =============================================================================
"""
if __name__ == '__main__':
    import langmuir_models2 as lm
    import data_slicer as ds
    
    print('processing shots')
    # from Workspace.shot_times import remaining_shots_stepan as shots_times
    from analysis_requests import ana_req_daihong as shots_times

    # shots_times = [('20180807.016', [0.0, 6.5])]
    
    model = lm.DoubleWeinlichSheathExpansion
#    fobj = open('C:/Users/luru/Desktop/GitHub/Data/log.txt', 'w')
    component_list = [ut.component_dict['UD'], ut.component_dict['LD']]
    probes = np.arange(1, 11)
    log = []
    for s, t in shots_times:
        # if s in done:
        # print(s+' already done!')
        # continue
        for component in component_list:
            try:
                print(f'Working on {s} {component}')
                log.append('---------------------\n')
                m = Manager(s, t, from_cache=False)
                m.load_probes(component,
                              probes,
                              force_supply_voltages=True,
                              load_processed=None)
                slicer_probe = m[component, 10]
                slicer = ds.VmaxClusterSlicer(slicer_probe, half_periods=2)
                slicer.remove_short_intervals(500)
                key = model, slicer
                # m.save()
                m.do_all('fit',
                         *key,
                         timed=True,
                         progressbar=True,
                         mute_warnings=[['ignore', r'invalid value'],
                                        ['ignore', r'divide by zero']],
                         verbose=2)
                m.saveStuff(key, component, quantities=['Te','ne'], directory=r'//sv-e4-fs-1/E4-Mitarbeiter/E4 Diagnostics/QRP_Langmuir_Sonden/OP1.2_TDU-Sonden/test analysis')
                # m.displayStuff(key,component,probes=probes)
                print(f'Done with {s}\n')
                log.append(f'Done with {s}\n')
                del m
            except Exception as e:
                print('---------------------')
                print(f'Error in {s}: {e}')
                log.append(f'Error in {s}: {e}\n')
                print('---------------------')
 #   for i in log:
 #       fobj.write(i)
 #   fobj.close()
"""