import subprocess
import os
import tempfile
import shutil
from pyrocko.snuffling import Snuffling, Param, Switch
from pyrocko import util, gui_util, guts
from xmlMarker import *
from PyQt4.QtCore import QUrl 
from PyQt4.QtGui import QDesktopServices

def get_magnitude(event):
    if event.magnitude:
        mag = event.magnitude
    elif event.moment_tensor:
        mag = event.moment_tensor.magnitude
    else:
        mag = 0.
    return mag

class MapMaker(Snuffling):
    '''
    <h1>Create a map containing event and station locations using googlemaps.</h1>

    Invokes the standard browser. Some browsers do not allow javascript to open
    and read the xml-file containing the necessary information due to the "Some-Origin-Policy".
    In that case you need to reset your standard browser.
    I.e.: Firefox on Linux do:
        <p><b>xdg-settings set default-web-browser firefox.desktop</b></p>
    <p>
    The plate boundary database is based on the work done by Peter Bird, who kindly permitted usage. <br>
    See: <i>50. Bird, P. (2003) An updated digital model of plate boundaries, Geochemistry Geophysics Geosystems, 4(3), 1027, doi:10.1029/2001GC000252.</i><br>

    Also available at <a href="http://peterbird.name/publications/2003_PB2002/2003_PB2002.htm">http://www.peterbird.name</a></p>
    <p>Please note, that in the current implementation the orogens (cross-hatched areas in <a href="http://peterbird.name/publications/2003_PB2002/Figure_01.gif">figure 1</a>) are not distinguished from plate boundaries.  The orogens are meant to mark areas where
    the plate model is known to be incomplete (and/or inapplicable).<br>
    This matter will be pointed out in future releases of this snuffling. 
    '''
    def setup(self):
        self.set_name('Create Map in GoogleMaps')
        self.add_parameter( Switch('Only Active Event', 'only_active', False))

    def call(self):
        self.cleanup()
        
        viewer = self.get_viewer()
        markers = viewer.get_markers()
        while True:
            try:
                active_event, active_stations = \
                        self.get_active_event_and_stations()
                break
            except AttributeError:
                if self.only_active == True:
                    self.set_parameter('only_active', False)
                    self.fail('No active event found.')

                active_stations = []
                active_event = None
                print 'Presumably no active event set'
                break

        station_list=[]
        stats = active_stations if self.only_active else viewer.stations.values() 

        for stat in stats:
            if not util.match_nslc(viewer.blacklist, stat.nsl()):
                xml_station_marker = XMLStationMarker(nsl = '.'.join(stat.nsl()),
                                                longitude = stat.lon,
                                                latitude = stat.lat,
                                                active = 'yes')

            station_list.append(xml_station_marker)
        active_station_list = StationMarkerList(stations=station_list)

        ev_marker_list = []
        if active_event is not None:

            depth = active_event.depth
            if depth is None:
                depth = 0.0

            xml_active_event_marker = XMLEventMarker(eventname = \
                    active_event.name,
                latitude=active_event.lat,
                longitude=active_event.lon,
                origintime=util.time_to_str(active_event.time),
                depth=depth,
                magnitude=get_magnitude(active_event),
                active='yes')
            ev_marker_list.append(xml_active_event_marker)    

        for m in markers:
            if self.only_active:
                break
            if isinstance(m, gui_util.EventMarker):
                ev = m.get_event()
                try:
                    depth = ev.depth
                    if depth is None:
                        depth = 0.0

                    xmleventmarker = XMLEventMarker(eventname=ev.name,
                                            longitude=ev.lon, 
                                            latitude=ev.lat, 
                                            origintime=util.time_to_str(ev.time), 
                                            depth=depth,
                                            magnitude=get_magnitude(ev),
                                            active='no')
                except guts.ValidationError:
                    print 'invalid format'
                    raise
            else:
                continue

            ev_marker_list.append(xmleventmarker)

        event_list=EventMarkerList(events= ev_marker_list)
        event_station_list = MarkerLists(station_marker_list=active_station_list,
                                         event_marker_list=event_list)
        event_station_list.validate()

        tempdir = tempfile.mkdtemp(prefix='googleMapsSnuffling-')

        for entry in ['loadxmldoc.js', 'map.html']:
            shutil.copy(os.path.join(self.module_dir(), entry), os.path.join(\
                    tempdir, entry))

        markers_fn = os.path.join(tempdir, 'markers.xml')
        dump_xml(event_station_list, filename=markers_fn)

        QDesktopServices.openUrl(QUrl('file://' + os.path.join(tempdir, \
                'map.html')))


def __snufflings__():
    return [ MapMaker() ]
