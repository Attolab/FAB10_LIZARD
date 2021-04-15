from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer, pyqtSlot, QThread
from pyqtgraph.dockarea import Dock
from pymodaq.daq_utils.daq_utils import linspace_step, PIDModelGeneric, \
    check_modules
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils import ellipses
from pymodaq.daq_utils.h5saver import H5Saver
import time
from datetime import datetime
import numpy as np
from collections import OrderedDict


class PIDModelLIZARD_FAB10(PIDModelGeneric):
    params = [
        {'title': 'Stabilization', 'name': 'stabilization', 'type': 'group',
         'children': [
            {'title': 'Ellipsis phase (rad):', 'name': 'ellipsis_phase',
             'type': 'float', 'value': 0},
            {'title': 'Error (rad):', 'name': 'error',
             'type': 'float', 'value': 0},
            {'title': 'Actuator absolute position (microns):',
             'name': 'actuator_absolute_position',
             'type': 'float', 'value': 0},
            {'title': 'Relative order to actuator (microns):',
             'name': 'relative_order_to_actuator',
             'type': 'float', 'value': 0},
            {'title': 'Laser wavelength (microns):', 'name': 'laser_wl',
             'type': 'float', 'value': 0.8},
            {'title': 'Correction sign:', 'name': 'correction_sign',
             'type': 'int', 'value': 1, 'min': -1, 'max': 1, 'step': 2},
            {'title': 'Offsets', 'name': 'offsets', 'type': 'group',
             'children': [
                {'title': 'Phase time zero (rad):', 'name': 'time_zero_phase',
                 'type': 'float', 'value': 0},
                {'title': 'Piezo (nm):', 'name': 'piezo_offset',
                 'type': 'float', 'value': 0}]},
         ]},
        {'title': 'Calibration', 'name': 'calibration', 'type': 'group',
         'expanded': True, 'visible': True, 'children': [
            {'title': 'Do calibration:', 'name': 'do_calibration',
             'type': 'bool', 'value': False},
            {'title': 'Timeout:', 'name': 'timeout', 'type': 'int',
             'value': 10000},
            {'title': 'Calibration', 'name': 'calibration_move',
             'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'title': 'Start pos:', 'name': 'start', 'type': 'float',
                 'value': -9544.},
                {'title': 'Stop pos:', 'name': 'stop', 'type': 'float',
                 'value': -9542.},
                {'title': 'Step size:', 'name': 'step', 'type': 'float',
                 'value': 0.04},
                {'title': 'Averaging:', 'name': 'average', 'type': 'int',
                 'value': 1}]},
            {'title': 'Ellipse params', 'name': 'calibration_ellipse',
             'type': 'group', 'expanded': True, 'visible': True, 'children': [
                 {'title': 'Dx:', 'name': 'dx', 'type': 'float', 'value': 0.0},
                 {'title': 'Dy:', 'name': 'dy', 'type': 'float', 'value': 0.0},
                 {'title': 'x0:', 'name': 'x0', 'type': 'float',
                  'value': 0.00},
                 {'title': 'y0:', 'name': 'y0', 'type': 'float',
                  'value': 0.00},
                 {'title': 'theta (°):', 'name': 'theta', 'type': 'float',
                  'value': 0.00}]}
         ]},
        {'title': 'Stabilized scan', 'name': 'stabilized_scan',
         'type': 'group', 'expanded': True, 'visible': True, 'children': [
            {'title': 'Do stabilized scan:', 'name': 'do_stabilized_scan',
             'type': 'bool', 'value': False},
            {'title': 'Timeout (ms):', 'name': 'timeout', 'type': 'int',
             'value': 5000},
            {'title': 'Stabilized scan parameters',
             'name': 'stabilized_scan_parameters', 'type': 'group',
             'expanded': True, 'visible': True, 'children': [
                {'title': 'Length (microns):', 'name': 'length',
                 'type': 'float', 'value': 1},
                {'title': 'Step size (microns):', 'name': 'step_size',
                 'type': 'float', 'value': 0.02},
                {'title': 'Iterations per setpoint:',
                 'name': 'iterations_per_setpoint', 'type': 'int',
                 'value': 2}]}
         ]}
    ]

    actuators = ['SmarActMCS']
    actuators_name = ['DELAY LINE']

    detectors_type = ['DAQ1D']
    detectors = ['LecroyWaverunner6Zi']
    detectors_name = ['OSCILLOSCOPE']

    check_modules(detectors, detectors_type, actuators)

    def __init__(self, pid_controller):
        super().__init__(pid_controller)

        self.dock_calib = Dock('Calibration')
        widget_calib = QtWidgets.QWidget()
        self.viewer_calib = Viewer1D(widget_calib)
        widget_ellipse = QtWidgets.QWidget()
        self.viewer_ellipse = Viewer1D(widget_ellipse)
        self.viewer_ellipse.show_data([np.zeros((10,)), np.zeros((10,))])
        self.dock_calib.addWidget(widget_calib)
        self.dock_calib.addWidget(widget_ellipse, row=0, col=1)
        self.pid_controller.dock_area.addDock(self.dock_calib)
        self.dock_calib.float()

        self.curr_time = time.perf_counter()
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        #  self.timer.timeout.connect(self.timeout)
        self.lsqe = ellipses.LSqEllipse()

        self.h5_saver = H5Saver()
        self.channel_pytables_array = None
        self.channel_pytables_array_time = None
        self.current_scan_path = None
        self.scope_spectrum_shape = None
        self.calibration_scan_shape = None
        self.stabilized_scan_shape = None

        self.calibration_scan_running = False
        self.stabilized_scan_running = False
        self.det_done_flag = False
        self.move_done_flag = False
        self.timeout_scan_flag = False
        self.calibration_done = False
        self.time_zero_phase_defined = False

        # The time zero corresponds to the temporal overlap of the XUV and IR
        # pulses. More precisely here it is defined as the starting position of
        # the delay stage during a calibration scan. Thus this position is
        # related to a particular calibration scan and defined by the user.
        # The time_zero_phase is set at the first iteration of the feedback
        # loop and will be offset to define the delay_phase which is then zero
        # at time zero.
        # The delay phase is unwrapped to keep track of the past from the
        # begining of the feedback loop. Thus there is a direct proportional
        # relation between delay_phase and the pump/probe delay.
        self.time_zero_phase = 0  # within [-pi, +pi]
        self.ellipsis_phase = 0  # within [-pi, +pi]
        self.delay_phase = 0  # unwrapped phase
        # This attribute is used to record every phase measurement (within
        # [-pi, +pi]). One can track the history of the evolution of the
        # phase using the unwrap function (see convert_input).
        self.phase_buffer = np.array([])

        self.actuator_absolute_position = 0
        self.setpoint = 0
        self.error = 0

    def ini_model(self):
        """ Initialize the PID model

        Defines all the action to be performed on the initialized modules
        (main pid, actuators, detectors).
        Either here for specific things (ROI, ...) or within the preset of the
        current model.
        """
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'output_limits',
            'output_limit_min').setValue(0)  # in microns
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'output_limits',
            'output_limit_min_enabled').setValue(False)
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'output_limits',
            'output_limit_max').setValue(15)  # in microns
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'output_limits',
            'output_limit_max_enabled').setValue(False)

        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'pid_constants',
            'kp').setValue(0.4)
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'pid_constants',
            'ki').setValue(0.05)
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'pid_constants',
            'kd').setValue(0)

        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'filter',
            'filter_enable').setValue(False)
        self.pid_controller.setpoint_sb.setValue(0)
        self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'sample_time').setValue(50)

        # Initialize a h5 file to save the data
        self.h5_saver.init_file(custom_naming=True)

    def get_ellipse_fit(self, center, width, height, theta):

        t = np.linspace(0, 2*np.pi, 1000)
        ellipse_x = (center[0]
                     + width*np.cos(t)*np.cos(theta)
                     - height*np.sin(t)*np.sin(theta))
        ellipse_y = (center[1]
                     + width*np.cos(t)*np.sin(theta)
                     + height*np.sin(t)*np.cos(theta))

        return ellipse_x, ellipse_y

    def update_settings(self, param):
        """Manage the change of a parameter on the UI.

        Get a parameter instance whose value has been modified by a user on the
        UI.

        Parameters
        ----------
        param: (Parameter) instance of Parameter object
        """
        if param.name() == 'do_calibration':
            if param.value():
                self.calibration_scan_running = True
                self.do_calibration()
                QtWidgets.QApplication.processEvents()
            else:
                self.calibration_scan_running = False
        if param.name() == 'do_stabilized_scan':
            if param.value():
                self.stabilized_scan_running = True
                self.do_stabilized_scan()
                QtWidgets.QApplication.processEvents()
            else:
                self.stabilized_scan_running = False

    def do_calibration(self):
        """ Perform a calibration scan

        This method is called by ticking the "do_calibration" parameter of the
            GUI.
        It should be called after the ROIs have been properly defined.
        It performs a scan of the delay line position and acquire the values of
            the modulated signals provided by the ROIs.
        The start position of this scan defines the time zero (the
            temporal superposition of the IR and XUV pulses).
        At the end of the scan the parameters of the ellipsis that fit the
            modulated signals in XY representation are set in the parameters
            and a window displays the result.
        The phase at time zero is also set.
        At the end of the scan the delay line returns to the time zero (the
            start position).
        """

        self.pid_controller.log_signal.emit(
            'A calibration scan has been launched !')

        # Create a new group "Calibration scan" is the h5 file
        self.h5_saver.add_scan_group(
            title='calibration_scan')
        self.current_scan_path = self.h5_saver.current_scan_group._v_pathname

        steps = linspace_step(self.settings.child('calibration',
                                                  'calibration_move',
                                                  'start').value(),
                              self.settings.child('calibration',
                                                  'calibration_move',
                                                  'stop').value(),
                              self.settings.child('calibration',
                                                  'calibration_move',
                                                  'step').value())

        self.calibration_scan_shape = [len(steps)]

        # Save the scan axis in the scan group
        if not self.h5_saver.is_node_in_group(self.h5_saver.current_scan_group,
                                              'scan_x_axis'):
            x_axis_meta = dict(
                units='micrometers',
                label='piezo device reading. x2 to get the delay between'
                      'pulses')
            self.h5_saver.add_navigation_axis(steps,
                                              self.h5_saver.current_scan_group,
                                              axis='x_axis',
                                              metadata=x_axis_meta)

        # Create a detector group for the oscilloscope
        self.h5_saver.add_det_group(
            where=self.h5_saver.current_scan_group,
            title='oscilloscope')
        detector_group = self.h5_saver.current_group

        # Create a data group 0D to save the ROIs
        self.h5_saver.add_data_group(
            where=detector_group,
            group_data_type='data0D',
            title='ROIs'
        )
        data_0d_group = self.h5_saver.current_group

        # Create a data group 1D to save the oscilloscope spectra
        self.h5_saver.add_data_group(
            where=detector_group,
            group_data_type='data1D',
            title='Oscilloscope'
        )
        data_1d_group = self.h5_saver.current_group

        # Create a channel group to store each ROI and the current time
        # ROI_00: first modulated signal
        # ROI_01: second modulated signal
        # ROI_02: signal offset
        # SECONDS: current time on the machine (seconds)
        # MINUTES: current time on the machine (minutes)
        # HOURS: current time on the machine (hours)
        # self.h5_saver.add_CH_group(
        #     where=data_0d_group,
        #     title='ROI_00'
        # )
        # self.h5_saver.add_CH_group(
        #     where=data_0d_group,
        #     title='ROI_01'
        # )
        # self.h5_saver.add_CH_group(
        #     where=data_0d_group,
        #     title='ROI_02'
        # )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='SECONDS'
        )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='MINUTES'
        )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='HOURS'
        )

        # Create a channel group to store the oscilloscope spectra in the
        # data1D group
        self.h5_saver.add_CH_group(
            where=data_1d_group,
            title='Time of flight (TOF) raw spectrum'
        )

        raw_spectrum_from_scope = []

        detector_data = np.zeros((steps.size, 2))
        current_time_array_seconds = np.zeros((steps.size, 1))
        current_time_array_minutes = np.zeros((steps.size, 1))
        current_time_array_hours = np.zeros((steps.size, 1))

        for mod in self.pid_controller.actuator_modules:
            mod.move_done_signal.connect(self.move_done)
        for mod in self.pid_controller.detector_modules:
            mod.grab_done_signal.connect(self.det_done)

        self.viewer_calib.x_axis = steps

        # The oscilloscope module (DAQ_Viewer object)
        scope_daq_viewer = self.pid_controller.detector_modules[0]
        # The actuator module (DAQ_Move object)
        delay_stage_daq_move = self.pid_controller.actuator_modules[0]

        # Use the ActiveDSO ClearDevice method.
        # This method will send a device clear signal to the oscilloscope.
        # Any unread response currently in the device output buffer will be
        # cleared. See Lecroy ActiveDSO Developer’s Guide.
        scope_daq_viewer.controller.DeviceClear(False)

        QtWidgets.QApplication.processEvents()
        QThread.msleep(1000)
        QtWidgets.QApplication.processEvents()

        for ind_step, step in enumerate(steps):
            self.move_done_flag = False
            delay_stage_daq_move.move_Abs(step)
            self.wait_for_move_done()

            self.det_done_flag = False
            scope_daq_viewer.grab_data()
            self.wait_for_det_done()

            raw_spectrum_from_scope = []
            raw_spectrum_from_scope = scope_daq_viewer.\
                data_to_save_export['data1D']['CH000']['data']

            # Initialize the PyTables array to save the data at the first
            # step of the scan.
            if ind_step == 0:
                # We store the shape of the data taken from the scope.
                # It will be used to initialize the saving of stabilized
                # scans
                self.scope_spectrum_shape = raw_spectrum_from_scope

                # Get the time axis of the scope
                scope_time_axis = scope_daq_viewer.\
                    data_to_save_export["data1D"]["CH000"]["x_axis"]["data"]

                self.channel_pytables_array = self.h5_saver.add_data(
                    channel_group=self.h5_saver.get_set_group(
                        where=self.current_scan_path+'/Detector000/Data1D',
                        name='Ch000'),
                    scan_type='scan1D',
                    scan_shape=self.calibration_scan_shape,
                    data_dict=dict(
                        data=self.scope_spectrum_shape,
                        x_axis=scope_time_axis),
                    title='TOF raw spectra',
                    enlargeable=False,
                    init=True,
                    add_scan_dim=True)

            self.channel_pytables_array.__setitem__(
                key=ind_step,
                value=raw_spectrum_from_scope)

            QtWidgets.QApplication.processEvents()
            QThread.msleep(300)
            QtWidgets.QApplication.processEvents()

            # The signal offset is taken from ROI_02 (mean value) which
            # should be out of any electron signal
            signal_offset = scope_daq_viewer.\
                data_to_save_export['data0D']['Measure_002']['data']
            spectrum_without_offset = []
            spectrum_without_offset = [elt - signal_offset for elt in
                                       raw_spectrum_from_scope]
            # The xuv_flux_normalization corresponds to the integral of the all
            # spectrum (i.e. the total number of photoelectrons).
            # Dividing by this value, we suppress the effect of the
            # fluctuations of the intensity of the XUV source.
            xuv_flux_normalization = np.trapz(spectrum_without_offset)

            # Measurement from ROI_00 (mean)
            # Correspond to the value of the first modulated signal m1
            # raw_m1 is the sum of the values in the ROI. Thus the ROI should
            # be defined properly by the user.
            raw_m1 = scope_daq_viewer.\
                data_to_save_export['data0D']['Measure_000']['data']
            # The width of the ROIs are needed. They are not equals since we
            # work with a time of flight signal.
            # region_m1 = scope_daq_viewer.ui.viewers[0].\
            #     roi_manager.ROIs["ROI_00"].pos()
            # width_m1 = region_m1[1] - region_m1[0]
            # We remove the offset times the width of the ROI to get a value
            # that is proportional to the number of electrons in this band
            # m1 = (raw_m1 - signal_offset*width_m1)/xuv_flux_normalization
            m1 = (raw_m1 - signal_offset) / xuv_flux_normalization

            # Measurement from ROI_01 (mean)
            raw_m2 = scope_daq_viewer.\
                data_to_save_export['data0D']['Measure_001']['data']
            # The width of the ROIs are needed. They are not equals since we
            # work with a time of flight signal.
            # region_m2 = scope_daq_viewer.ui.viewers[0].\
            #     roi_manager.ROIs["ROI_01"].pos()
            # width_m2 = region_m2[1] - region_m2[0]
            # m2 = (raw_m2 - signal_offset*width_m2)/xuv_flux_normalization
            m2 = (raw_m2 - signal_offset) / xuv_flux_normalization

            detector_data[ind_step, 0] = m1
            detector_data[ind_step, 1] = m2

            current_time = datetime.now().time()
            current_time_array_seconds[ind_step] = current_time.second
            current_time_array_minutes[ind_step] = current_time.minute
            current_time_array_hours[ind_step] = current_time.hour

            if not self.calibration_scan_running:
                break

        self.viewer_calib.show_data([detector_data[:, 0],
                                     detector_data[:, 1]])

        self.lsqe.fit([detector_data[:, 0],
                       detector_data[:, 1]])
        center, width, height, theta = self.lsqe.parameters()
        ellipse_x, ellipse_y = self.get_ellipse_fit(center, width, height,
                                                    theta)

        self.viewer_ellipse.plot_channels[0].\
            setData(x=detector_data[:, 0],
                    y=detector_data[:, 1])
        self.viewer_ellipse.plot_channels[1].setData(x=ellipse_x,
                                                     y=ellipse_y)

        self.settings.child('calibration', 'calibration_ellipse', 'x0').\
            setValue(center[0])
        self.settings.child('calibration', 'calibration_ellipse', 'y0').\
            setValue(center[1])
        self.settings.child('calibration', 'calibration_ellipse', 'dx').\
            setValue(width)
        self.settings.child('calibration', 'calibration_ellipse', 'dy').\
            setValue(height)
        self.settings.child('calibration', 'calibration_ellipse',
                            'theta').setValue(np.rad2deg(theta))

        QtWidgets.QApplication.processEvents()

        # Save the data of the two modulated signals
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch000'),
            data_dict=dict(data=detector_data[:, 0]),
            title='Modulated signal m1 (averaged, offset, normalized)')
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch001'),
            data_dict=dict(data=detector_data[:, 1]),
            title='Modulated signal m2 (averaged, offset, normalized)')

        # Save the current time (seconds)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch003'),
            data_dict=dict(data=current_time_array_seconds),
            title='Current time (seconds)')
        # Save the current time (minutes)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch004'),
            data_dict=dict(data=current_time_array_minutes),
            title='Current time (minutes)')
        # Save the current time (hours)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch005'),
            data_dict=dict(data=current_time_array_hours),
            title='Current time (hours)')

        # The actuator should go back to the starting position (time zero) at
        # the end of the calibration scan.
        self.move_done_flag = False
        delay_stage_daq_move.move_Abs(steps[0])
        self.wait_for_move_done()

        self.calibration_done = True
        self.settings.child('calibration', 'do_calibration').setValue(False)

        for mod in self.pid_controller.actuator_modules:
            mod.move_done_signal.disconnect(self.move_done)
        for mod in self.pid_controller.detector_modules:
            mod.grab_done_signal.disconnect(self.det_done)

        self.pid_controller.log_signal.emit(
            'The calibration scan is finished! '
            'You can INIT, PLAY and uncheck PAUSE to launch the feedback loop')

    def do_stabilized_scan(self):
        """Perform a stabilized scan.

        It will start at the current position (to prevent from too big change
            of setpoint value).
        By default it will increase the setpoint.
        The user should give the length, the step size and the number of PID
            loop iterations on each step.
        This method should be called while the PID loop is running.
        """
        # Check that the calibration has been done and that the PID loop is
        # running. If not the function is not executed.
        if (not self.calibration_done
                or self.pid_controller.PIDThread.pid_runner.paused):
            self.pid_controller.log_signal.emit(
                'The PID loop should run before launching a stabilized scan! '
                'Push the PAUSE button!')
            return

        self.pid_controller.log_signal.emit(
            'A stabilized scan has been launched !')

        # Number of PID loop iterations on each step
        iterations_per_setpoint = self.settings.child(
            'stabilized_scan',
            'stabilized_scan_parameters',
            'iterations_per_setpoint').value()
        # Scan length in microns
        scan_length = self.settings.child('stabilized_scan',
                                          'stabilized_scan_parameters',
                                          'length').value()
        # Scan step in microns
        scan_step = self.settings.child('stabilized_scan',
                                        'stabilized_scan_parameters',
                                        'step_size').value()
        # Laser wavelength in microns
        laser_wl = self.settings.child('stabilization', 'laser_wl').value()

        setpoint_steps = linspace_step(
            self.delay_phase,
            self.delay_phase + 8*np.pi*scan_length/laser_wl,
            8*np.pi*scan_step/laser_wl
        )

        # The shape of the scan will be used to initialize the PyTables array
        # in convert_input method
        self.stabilized_scan_shape = [
            len(setpoint_steps)*iterations_per_setpoint]

        stabilized_scan_metadata = dict(
            iterations_per_setpoint=str(iterations_per_setpoint)
        )
        self.h5_saver.add_scan_group(
            title='stabilized scan',
            metadata=stabilized_scan_metadata)
        self.current_scan_path = self.h5_saver.current_scan_group._v_pathname

        # Save the setpoint steps of the scan
        setpoint_steps_metadata = dict(units='radians')
        self.h5_saver.add_navigation_axis(setpoint_steps,
                                          self.h5_saver.current_scan_group,
                                          axis='x_axis',
                                          metadata=setpoint_steps_metadata)

        # Create a detector group for the oscilloscope
        self.h5_saver.add_det_group(
            where=self.h5_saver.current_scan_group,
            title='oscilloscope')
        detector_group = self.h5_saver.current_group

        # Create a data group 0D to save the ROIs
        self.h5_saver.add_data_group(
            where=detector_group,
            group_data_type='data0D',
            title='current time'
        )
        data_0d_group = self.h5_saver.current_group

        # Create a channel group to store the current time
        # SECONDS: current time on the machine (seconds)
        # MINUTES: current time on the machine (minutes)
        # HOURS: current time on the machine (hours)
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='SECONDS'
        )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='MINUTES'
        )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='HOURS'
        )

        # Create channel groups to store the delay phase and the error at each
        # iteration of the feedback loop
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='DELAY PHASE'
        )
        self.h5_saver.add_CH_group(
            where=data_0d_group,
            title='ERROR'
        )

        # Create a data group 1D to save the oscilloscope spectra
        self.h5_saver.add_data_group(
            where=detector_group,
            group_data_type='data1D',
            title='Oscilloscope spectra'
        )
        data_1d_group = self.h5_saver.current_group

        # Create a channel group to store the oscilloscope spectra in the
        # data1D group
        self.h5_saver.add_CH_group(
            where=data_1d_group,
            title='Time of flight (TOF) raw spectrum'
        )

        raw_spectrum_from_scope = []

        for mod in self.pid_controller.detector_modules:
            mod.grab_done_signal.connect(self.det_done)

        # Initialize the arrays to save the clock of the computer at each
        # iteration
        array_size = setpoint_steps.size*iterations_per_setpoint
        current_time_array_seconds = np.zeros((array_size, 1))
        current_time_array_minutes = np.zeros((array_size, 1))
        current_time_array_hours = np.zeros((array_size, 1))

        # Initialize the arrays to store the delay phase and the error at each
        # iteration
        delay_phase_array = np.zeros((array_size, 1))
        error_array = np.zeros((array_size, 1))

        for ind_setpoint, setpoint in enumerate(setpoint_steps):
            self.pid_controller.setpoint_sb.setValue(setpoint)

            for iteration in range(iterations_per_setpoint):
                self.det_done_flag = False
                self.wait_for_det_done()

                raw_spectrum_from_scope = []
                raw_spectrum_from_scope = self.pid_controller.\
                    detector_modules[0].\
                    data_to_save_export['data1D']['CH000']['data']

                # Initialize the PyTables array to save the data at the first
                # step of the scan.
                if ind_setpoint == 0 and iteration == 0:
                    # We store the shape of the data taken from the scope.
                    # It will be used to initialize the saving of stabilized
                    # scans
                    self.scope_spectrum_shape = raw_spectrum_from_scope

                    # Get the time axis of the scope
                    scope_time_axis = self.pid_controller.detector_modules[0]. \
                        data_to_save_export["data1D"]["CH000"]["x_axis"]["data"]

                    # Initialization of the pytable array for the 1D data of
                    # the scope
                    self.channel_pytables_array = self.h5_saver.add_data(
                        channel_group=self.h5_saver.get_set_group(
                            where=self.current_scan_path+'/Detector000/Data1D',
                            name='Ch000'),
                        scan_type='scan1D',
                        scan_shape=self.stabilized_scan_shape,
                        data_dict=dict(
                            data=self.scope_spectrum_shape,
                            x_axis=scope_time_axis),
                        title='TOF raw spectra',
                        enlargeable=False,
                        init=True,
                        add_scan_dim=True)

                current_iteration = (ind_setpoint*iterations_per_setpoint
                                     + iteration)

                self.channel_pytables_array.__setitem__(
                    key=current_iteration,
                    value=raw_spectrum_from_scope)

                # Save the current time
                current_time = datetime.now().time()
                current_time_array_seconds[current_iteration] = current_time.second
                current_time_array_minutes[current_iteration] = current_time.minute
                current_time_array_hours[current_iteration] = current_time.hour

                # Save the current delay phase and error
                delay_phase_array[current_iteration] = self.delay_phase
                error_array[current_iteration] = self.error

                QtWidgets.QApplication.processEvents()
                QThread.msleep(300)
                QtWidgets.QApplication.processEvents()

            if not self.stabilized_scan_running:
                break

        # Save the current time (seconds)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch003'),
            data_dict=dict(data=current_time_array_seconds),
            title='Current time (seconds)')

        # Save the current time (minutes)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch004'),
            data_dict=dict(data=current_time_array_minutes),
            title='Current time (minutes)')

        # Save the current time (hours)
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch005'),
            data_dict=dict(data=current_time_array_hours),
            title='Current time (hours)')

        # Save the delay phase
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch006'),
            data_dict=dict(data=delay_phase_array),
            title='Delay phase (radians)')

        # Save the error
        self.h5_saver.add_data(
            channel_group=self.h5_saver.get_set_group(
                where=self.current_scan_path+'/Detector000/Data0D',
                name='Ch007'),
            data_dict=dict(data=error_array),
            title='Error (radians)')

        self.stabilized_scan_running = False
        # Setting calibration_done to False, we force the user to redo a fresh
        # calibration scan before launching a stabilized scan.
        self.calibration_done = False

        for mod in self.pid_controller.detector_modules:
            mod.grab_done_signal.disconnect(self.det_done)

        self.pid_controller.log_signal.emit(
            'The stabilized scan is finished! If you want to restart a '
            'new one, you have to do a fresh calibration scan.')

    def timeout(self):
        """
            Send the status signal *'Time out during acquisition'* and stop the
            timer.
        """
        self.timeout_scan_flag = True
        self.timer.stop()
        self.pid_controller.update_status("Timeout during acquisition",
                                          log_type='log')

    def wait_for_det_done(self):
        """Wait for all the detectors to be ready or timeout.
        """
        self.timeout_scan_flag = False
        self.timer.start(self.settings.child('calibration', 'timeout').value())
        while not(self.det_done_flag or self.timeout_scan_flag):
            # wait for grab done signals to end
            QtWidgets.QApplication.processEvents()
        self.timer.stop()

    pyqtSlot(OrderedDict)
    def det_done(self, data):
        """Called each time the oscilloscope finished his acquisition.

        Parameters
        ----------
        data: ??
        """

        self.det_done_flag = True

    def wait_for_move_done(self):
        """Wait for all the actuators to be ready or timeout.
        """
        self.timeout_scan_flag = False
        self.timer.start(self.settings.child('calibration', 'timeout').value())
        while not(self.move_done_flag or self.timeout_scan_flag):
            # wait for move done signals to end
            QtWidgets.QApplication.processEvents()
        self.timer.stop()

    pyqtSlot(str, float)
    def move_done(self, actuator_title, position):
        """Triggered by a signal from the actuator (DAQ_Move object)

        :param actuator_title: (str) title of the actuator module
        :param position: (float) position of the actuator in microns
        """
        self.actuator_absolute_position = position
        self.move_done_flag = True

    def get_phi_from_xy(self, x, y):
        """Return the angle corresponding to the point (x,y) on the ellipsis.

        :param x: (float) abscissa in a generic frame of reference. We choose
            it to be the first modulated signal M1.
        :param y: (float) ordinate in a generic frame of reference. We choose
            it to be the second modulated signal M2.
        :return: (float) corresponding angle on the ellipsis whose parameters
            has been set during the calibration scan (in radians).
        """
        x0 = self.settings.child('calibration', 'calibration_ellipse',
                                 'x0').value()
        y0 = self.settings.child('calibration', 'calibration_ellipse',
                                 'y0').value()
        dx = self.settings.child('calibration', 'calibration_ellipse',
                                 'dx').value()
        dy = self.settings.child('calibration', 'calibration_ellipse',
                                 'dy').value()
        theta = np.deg2rad(self.settings.child('calibration',
                                               'calibration_ellipse',
                                               'theta').value())

        # Apply a rotation of -theta to get the ellipsis on main axes
        v = np.array([x, y])
        rot = np.array(
            [[np.cos(-theta), -np.sin(-theta)],
             [np.sin(-theta), np.cos(-theta)]])
        vprim = rot @ v

        phi = np.arctan2((vprim[1]
                          + x0*np.sin(theta)
                          - y0*np.cos(theta))/dy,
                         (vprim[0]
                          - x0*np.cos(theta)
                          - y0*np.sin(theta))/dx)

        self.settings.child('stabilization', 'ellipsis_phase').setValue(phi)
        self.ellipsis_phase = phi

        return phi

    def convert_input(self, measurements):
        """Return a measured phase (delay) from the oscilloscope spectrum.

        Convert the measurements from the ROIs in the electronic spectrum to a
        measured phase in rad (same dimensionality as the setpoint). The output
        feeds the PID module.

        Parameters
        ----------
        measurements: (list) List of object from which the model extract a
        value of the same units as the setpoint

        Returns
        -------
        float: unwrapped phase in radians. The origin of the phase axis is the
                phase at time zero.
        """

        # Spectrum from the oscilloscope
        raw_spectrum_from_scope =\
            measurements['OSCILLOSCOPE']['data1D']['CH000']['data']

        # The signal offset is taken from ROI_02 (mean value), which should be
        # out of any electron signal.
        signal_offset =\
            measurements['OSCILLOSCOPE']['data0D']['Measure_002']['data']
        spectrum_without_offset =\
            [elt - signal_offset for elt in raw_spectrum_from_scope]
        # The modulated signals are normalized by the total number of
        # photoelectrons in the signal to prevent from XUV intensity
        # fluctuations.
        normalization_constant = np.trapz(spectrum_without_offset)

        # Measurement from ROI_00 (mean value)
        # Correspond to the value of the first modulated signal m1
        raw_m1 = measurements['OSCILLOSCOPE']['data0D']['Measure_000']['data']
        m1 = (raw_m1 - signal_offset)/normalization_constant
        # Measurement from ROI_01 (mean value)
        # Correspond to the value of the second modulated signal m2
        raw_m2 = measurements['OSCILLOSCOPE']['data0D']['Measure_001']['data']
        m2 = (raw_m2 - signal_offset)/normalization_constant

        # The phase returned by get_phi_from_xy is within [-pi, +pi].
        phi = self.get_phi_from_xy(m1, m2)

        # The first call of convert_input method is done after the user push
        # the PLAY button (launch of the PID loop).
        # It should be done just after the calibration scan. Which means that
        # we are at time zero.
        # That is how we define the phase at time zero.
        if not self.time_zero_phase_defined:
            self.time_zero_phase = phi
            self.settings.child('stabilization', 'offsets',
                                'time_zero_phase').setValue(
                self.time_zero_phase)
            self.time_zero_phase_defined = True

        # Record the measured phases. The phases in self.phase_buffer should be
        # within [-pi,+pi].
        self.phase_buffer = np.append(self.phase_buffer, [phi])

        # Perform the unwrap operation and offset by time_zero_phase. Thus
        # delay_phase is zero at time zero.
        phase_buffer_unwrapped = np.unwrap(
            self.phase_buffer - self.time_zero_phase)
        unwrapped_phase = phase_buffer_unwrapped[-1]
        self.delay_phase = unwrapped_phase

        return unwrapped_phase

    def convert_output(self, phase_correction, dt, stab=True):
        """
        Convert the phase correction from the PID module to an order in
        absolute value for the actuator.

        Parameters
        ----------
        phase_correction: (float) output value from the PID module.
                This phase correction in radians should be within [-pi,+pi].
        dt: ??
        stab: ??

        Returns
        -------
        list: the converted output as a list (if there are a few actuators).
                Absolute value in microns for the piezo actuator.
        """

        # Laser wavelength in microns
        laser_wl = self.settings.child('stabilization', 'laser_wl').value()
        # This parameter is used to easily change the sign of the correction
        correction_sign = self.settings.child('stabilization',
                                              'correction_sign').value()

        # Get the current position of the actuator
        self.actuator_absolute_position =\
            self.pid_controller.actuator_modules[0].current_position

        absolute_order_to_actuator = (
                self.actuator_absolute_position
                + correction_sign*phase_correction*laser_wl/(8*np.pi))

        self.settings.child(
            'stabilization',
            'actuator_absolute_position').setValue(
            self.actuator_absolute_position)

        # Get the value of the error (in rad)
        self.setpoint = self.pid_controller.settings.child(
            'main_settings', 'pid_controls', 'set_point').value()
        self.error = self.delay_phase - self.setpoint
        self.settings.child('stabilization', 'error').setValue(self.error)

        self.settings.child(
            'stabilization',
            'relative_order_to_actuator').setValue(
            absolute_order_to_actuator - self.actuator_absolute_position)

        return [absolute_order_to_actuator]


if __name__ == '__main__':
    pass
