#-*- coding: utf-8 -*-

import os

from pyrocko.snuffling import Choice, Snuffling, Switch
from pyrocko.pile_viewer import Marker

# TO-DO - save to VELEST format.

# Defaults to user's home folder.
USER_DEF_ROOT_FOLDER = \
os.path.expanduser('~')

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
        # self.add_parameter(Switch('Save in Velest format',
        #     'velest_format', False))

        self.setup_gui(reloaded=True)
        self.set_live_update(False)
    
    def parse_markers(self):
        self.cleanup()
        
        if self.save_asc_markers and not self.save_sel_markers:
            active_event, asc_phase_markers = \
                self.get_active_event_and_phase_markers()
            asc_phase_markers.insert(0, active_event)
            self.selected_markers = asc_phase_markers
            self.ae_name = active_event.get_event().name
        elif not self.save_asc_markers and self.save_sel_markers:
            self.selected_markers = self.get_selected_markers
        else:
            self.selected_markers = self.get_markers()

    def save_markers(self):
        if self.save_asc_markers and not self.save_sel_markers:
            if self.name_to_path and self.name_as_filename:
                save_path = self.root_folder + "/" + self.ae_name
                try:
                    os.makedirs(save_path)
                except OSError:
                    if not os.path.isdir(save_path):
                        raise
                self.save_path = save_path + "/" + self.ae_name
            elif self.name_to_path:
                save_path = self.root_folder + "/" + self.ae_name
                try:
                    os.makedirs(save_path)
                except OSError:
                    if not os.path.isdir(save_path):
                        raise
                self.save_path = save_path + "/picks"
            elif self.name_as_filename:
                self.save_path = self.root_folder + "/" + self.ae_name
            else:
                self.save_path = self.root_folder + "/picks"
        else:
            self.save_path = self.root_folder + "/picks"
        
        Marker.save_markers(self.selected_markers, self.save_path)
    
    def save_velest(self):
        pass
    
    def call(self):
        self.parse_markers()
        self.save_markers()
        print("Markers saved to: {}.".format(self.save_path))

def __snufflings__():
    """Returns a list of snufflings to be exported by this module."""

    return[QuickSave()]

