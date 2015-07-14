import subprocess
import os
import tempfile
import shutil
from pyrocko.snuffling import Snuffling, Param, Switch, Choice, NoViewerSet
from pyrocko import util, gui_util, guts, model
from xmlMarker import (XMLEventMarker, EventMarkerList, XMLStationMarker,
    StationMarkerList, MarkerLists, dump_xml)
from PyQt4.QtCore import QUrl 
from PyQt4.QtGui import QDesktopServices

g_counter = 0

def get_magnitude(event):
    if event.magnitude:
        mag = event.magnitude
    elif event.moment_tensor:
        mag = event.moment_tensor.moment_magnitude()
    else:
        mag = 0.
    return float(mag)

def convert_event_marker(marker):
    ev = marker.get_event()
    depth = ev.depth
    if depth is None:
        depth = 0.0
    ev_name = ev.name if ev.name else '(Event)'
    xmleventmarker = XMLEventMarker(eventname=ev_name,
                            longitude=float(ev.lon),
                            latitude=float(ev.lat),
                            origintime=util.time_to_str(ev.time),
                            depth=float(depth),
                            magnitude=float(get_magnitude(ev)),
                            active=['no', 'yes'][marker._active])

    return xmleventmarker

class MapMaker(Snuffling):
    '''
    <html>
    <body>
    <h1>Map event and stations with OpenStreetMap or Google Maps</h1>

    <p>
    Invokes the standard browser if "Open in external browser" is selected.
    Some browsers do not allow javascript to open and read the xml-file 
    containing the necessary information due to the "Same-Origin-Policy".
    In that case you need to reset your standard browser. I.e.: Firefox on
    Linux do: <tt>xdg-settings set default-web-browser firefox.desktop</tt>
    </p>
    <p>
    Clicking one of the plate boundary lines shows a reference regarding that
    plate boundary.
    </p>
    <p>
    The plate boundary database is based on the work done by Peter Bird, who kindly permitted usage. <br>
    See: <i>50. Bird, P. (2003) An updated digital model of plate boundaries, Geochemistry Geophysics Geosystems, 4(3), 1027, doi:10.1029/2001GC000252.</i>
    <br>
    Also available at 
    <a href="http://peterbird.name/publications/2003_PB2002/2003_PB2002.htm">http://www.peterbird.name</a>
    <br>
    Please note, that in the current implementation the orogens (cross-hatched
    areas in 
    <a href="http://peterbird.name/publications/2003_PB2002/Figure_01.gif">figure 1</a>)
    are not distinguished from plate boundaries.  The orogens are meant to 
    mark areas where the plate model is known to be incomplete (and/or inapplicable).<br>
    This matter will be pointed out in future releases of this snuffling. 
    </p>
    <p>
    This snuffling can also be called from the command line, if it is stored in
    the default pyrocko location under $HOME/.snufflings<br>
    e.g.:
    <code>
python $HOME/.snufflings/map/snuffling.py --stations=stations.pf
--events=events_test.pf
</code>
</p>
    </body>
    </html>
    '''
    def setup(self):
        self.set_name('Map')
        self.add_parameter(Switch('Only active event', 'only_active', False))
        self.add_parameter(Switch('Open in external browser',
                                  'open_external', False))
        self.add_parameter(Choice('Provider', 'map_kind', 'OpenStreetMap',
                                  ['OpenStreetMap', 'Google Maps']))

        self.set_live_update(False)

    def call(self):
        self.cleanup()

        try:
            viewer = self.get_viewer()
            cli_mode = False
        except NoViewerSet:
            viewer = None
            cli_mode = True
        
        if not cli_mode:
            if self.only_active:
                active_event, active_stations = \
                self.get_active_event_and_stations()
            else:
                active_event = None
                active_stations = viewer.stations.values()
        elif cli_mode:
            active_stations = self.stations

        station_list=[]
        if active_stations:
            for stat in active_stations:
                if (viewer and not util.match_nslc(viewer.blacklist, stat.nsl())) or cli_mode:
                    xml_station_marker = XMLStationMarker(
                        nsl='.'.join(stat.nsl()),
                        longitude = float(stat.lon),
                        latitude = float(stat.lat),
                        active = 'yes')

                    station_list.append(xml_station_marker)

        else:
            stations_list = []
        active_station_list = StationMarkerList(stations=station_list)

        if self.only_active:
            markers = [viewer.get_active_event_marker()]
        else:
            if not cli_mode:
                tmin, tmax = self.get_selected_time_range(fallback=True)
                markers = [m for m in viewer.get_markers()
                           if isinstance(m, gui_util.EventMarker) and\
                          m.tmin>=tmin and m.tmax<=tmax]

            else:
                markers = self.markers

        ev_marker_list = []
        for m in markers:
            xmleventmarker = convert_event_marker(m)
            ev_marker_list.append(xmleventmarker)

        event_list=EventMarkerList(events= ev_marker_list)
        event_station_list = MarkerLists(
            station_marker_list=active_station_list,
            event_marker_list=event_list)

        event_station_list.validate()

        tempdir = tempfile.mkdtemp(dir=self.tempdir())

        if self.map_kind == 'Google Maps':
            map_fn = 'map_googlemaps.html'
        elif self.map_kind == 'OpenStreetMap':
            map_fn = 'map_osm.html'

        url = 'file://' + tempdir + '/' + map_fn

        for entry in ['loadxmldoc.js', map_fn]:
            if cli_mode:
                snuffling_dir = os.environ['HOME']+'/.snufflings/map/'
            else:
                snuffling_dir = self.module_dir()

            shutil.copy(os.path.join(snuffling_dir, entry),
                    os.path.join(tempdir, entry))

        markers_fn = os.path.join(tempdir, 'markers.xml')
        dump_xml(event_station_list, filename=markers_fn)

        if self.open_external:
            QDesktopServices.openUrl(QUrl(url))
        else:
            global g_counter
            g_counter += 1
            self.web_frame(
                url,
                name='Map %i (%s)' % (g_counter, self.map_kind))


    def configure_cli_parser(self, parser):

         parser.add_option(
            '--events',
            dest='events_filename',
            default=None, 
            metavar='FILENAME',
            help='Read markers from FILENAME')

         parser.add_option(
            '--markers',
            dest='markers_filename',
            default=None, 
            metavar='FILENAME',
            help='Read markers from FILENAME')

         parser.add_option(
            '--stations',
            dest='stations_filename',
            default=None,
            metavar='FILENAME',
            help='Read stations from FILENAME')

         parser.add_option(
            '--provider',
            dest='map_provider',
            default='google',
            help='map provider [google | osm] (default=osm)')


def __snufflings__():
    return [ MapMaker() ]


if __name__ == '__main__':
    util.setup_logging('map.py', 'info')
    s = MapMaker()
    options, args, parser = s.setup_cli()
    s.markers = [] 

    if options.stations_filename:
        stations = model.load_stations(options.stations_filename)
        s.stations = stations
    else:
        s.stations = None

    if options.events_filename:
        events = model.load_events(filename=options.events_filename)
        markers = [gui_util.EventMarker(e) for e in events]
        s.markers.extend(markers)

    if options.markers_filename:
        markers = gui_util.load_markers(options.markers_filename)
        s.markers.extend(markers)
    s.open_external = True
    mapmap = {'google': 'Google Maps', 'osm': 'OpenStreetMap'}
    s.map_kind = mapmap[options.map_provider]
    s.call()
