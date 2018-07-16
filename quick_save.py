#-*- coding: utf-8 -*-

import datetime
import os

from pyrocko.snuffling import Choice, Snuffling, Switch
from pyrocko.pile_viewer import Marker

# Defaults to user's home folder. os.path.expanduser('~')
USER_DEF_ROOT_FOLDER = \
os.path.expanduser('~')


def folder_check(save_path):
    try:
        os.makedirs(save_path)
    except OSError:
        if not os.path.isdir(save_path):
            raise


class QuickSave(Snuffling):

    """

    <html>
    <h2>Quick Save</h2>
    <body>
    <p>
    Choose the predefined or user defined path from the <b>Root folder</b> 
    menu, desired <b>options</b> and press <b>Run</b>. Path to marker file is
    printed to the terminal.<br>
    User defined path is selected by default and can be modified in the source
    file. It is given as a string under <i>USER_DEF_ROOT_FOLDER</i> variable
    (line 11). Its default value points to user's home folder.
    </p>
    <p>
    Options:</br>
    <ul>
    <li>Save only selected markers.</li>
    <li>Save only markers associated to the active event.</li>
        <ul>
        <li>Add active event name to path. Creates new folder if needed.</li>
        <li>Use active event name as filename.</li>
        <li>Save to Velest format.</li>
        </ul>
    </ul>
    </p>
    <p>
    Author: G. Rajh.
    </p>
    </body>
    </html>

    """

    def setup(self):
        self.set_name('Quick Save')

        home_folder = os.path.expanduser('~')

        self.add_parameter(Choice('Root folder', 'root_folder',
            USER_DEF_ROOT_FOLDER, [USER_DEF_ROOT_FOLDER,
            home_folder + '/Desktop',
            home_folder + '/Documents'
            ]))

        self.add_parameter(Switch(
            'Save only selected markers', 'save_sel_markers', False))
        self.add_parameter(Switch(
            'Save only markers associated\nto the active event',
            'save_asc_markers', False))
        self.add_parameter(Switch('Add active event name to path',
            'name_to_path', False))
        self.add_parameter(Switch('Use active event name as filename',
            'name_as_filename', False))
        self.add_parameter(Switch('Save to Velest format',
            'velest_format', False))

        self.setup_gui(reloaded=True)
        self.set_live_update(False)
    
    def parse_markers(self):
        self.cleanup()
        
        if self.save_asc_markers and not self.save_sel_markers:
            active_event_marker, asc_phase_markers = \
                self.get_active_event_and_phase_markers()
            self.phase_markers = asc_phase_markers[:]
            asc_phase_markers.insert(0, active_event_marker)
            self.selected_markers = asc_phase_markers
            self.active_event = active_event_marker.get_event()
            self.ae_name = self.active_event.name

        elif not self.save_asc_markers and self.save_sel_markers:
            self.selected_markers = self.get_selected_markers()

        else:
            self.selected_markers = self.get_markers()

    def save_markers(self):
        if self.save_asc_markers and not self.save_sel_markers:

            if self.name_to_path and self.name_as_filename:
                save_path = self.root_folder + "/" + self.ae_name
                folder_check(save_path)
                self.save_path = save_path + "/" + self.ae_name

            elif self.name_to_path:
                save_path = self.root_folder + "/" + self.ae_name
                folder_check(save_path)
                self.save_path = save_path + "/picks"

            elif self.name_as_filename:
                self.save_path = self.root_folder + "/" + self.ae_name

            else:
                self.save_path = self.root_folder + "/picks"

        else:
            self.save_path = self.root_folder + "/picks"
        
        Marker.save_markers(self.selected_markers, self.save_path)
    
    def save_velest(self):
        start_time = datetime.datetime(1970, 1, 1)
        event_time = self.active_event.time
        event_date = start_time.utcfromtimestamp(event_time)
        event_year = int(str(event_date.year)[2:])
        event_month = event_date.month
        event_day = event_date.day
        event_hour = event_date.hour
        event_min = event_date.minute
        event_sec = event_date.second + (
            round(event_date.microsecond / 10**6, 2))
        event_lat = round(self.active_event.lat, 4)
        event_lon = round(self.active_event.lon, 4)
        event_depth = round(self.active_event.depth / 10**3, 2)
        event_mag = round(self.active_event.magnitude, 2)

        if event_lat >= 0:
            event_NS = "N"
        else:
            event_NS = "S"

        if event_lon >= 0:
            event_EW = "E"
        else:
            event_EW = "W"

        event_line = "{:02d}{:02d}{:02d} {:02d}{:02d} {:05.2f} {:7.4f}{} \
{:8.4f}{}{:8.2f}{:7.2f}     99  0.0 0.00  1.0  1.0 \n".format(
            event_year, event_month, event_day, event_hour, event_min,
            event_sec, event_lat, event_NS, event_lon, event_EW, event_depth,
            event_mag 
            )

        self.save_path_velest = self.save_path + "_velest"

        if self.name_to_path or self.name_as_filename:
            velest_file = open(self.save_path_velest, 'w')
        else:
            velest_file = open(self.save_path_velest, 'a+')

        velest_file.write(event_line)

        phase_i = 0
        phase_num = len(self.phase_markers)

        for phase_marker in self.phase_markers:
            phase_i += 1            
            phase_station = list(phase_marker.nslc_ids)[0][1]
            phase_name = phase_marker._phasename
            phase_tmin = phase_marker.tmin
            phase_tmax = phase_marker.tmax
            phase_unc = phase_marker._uncertainty

            if phase_tmin == phase_tmax:
                phase_time = phase_tmin
            else:
                phase_mid_int = (phase_tmax - phase_tmin) * 0.5
                phase_time = phase_tmin + phase_mid_int

            phase_rtime = round(phase_time - event_time, 2)

            if phase_unc <= 0.1:
                phase_unc_class = 0
            elif 0.1 < phase_unc <= 0.2:
                phase_unc_class = 1
            elif 0.2 < phase_unc <= 0.5:
                phase_unc_class = 2
            else:
                phase_unc_class = 3

            if phase_i % 6.0 == 0 or phase_i == phase_num:
                phase_block = "{:4s}{}{}{:06.2f}\n".format(phase_station,
                    phase_name, phase_unc_class, phase_rtime)
            else:
                phase_block = "{:4s}{}{}{:06.2f}".format(phase_station,
                    phase_name, phase_unc_class, phase_rtime)
                

            velest_file.write(phase_block)

        velest_file.write('\n')
        velest_file.close()
    
    def call(self):
        self.parse_markers()
        self.save_markers()
        print('Markers saved to: {}.'.format(self.save_path))

        if self.velest_format:
            self.save_asc_markers = True
            self.save_sel_markers = False
            self.parse_markers()
            self.save_velest()
            print('Velest format saved to: {}.'.format(self.save_path_velest))
        else:
            pass

def __snufflings__():
    """Returns a list of snufflings to be exported by this module."""

    return[ QuickSave() ]

