from pyrocko.snuffling import Snuffling, Param, Switch
from pyrocko import util, gui_util
from xmlMarker import *
from PyQt4.QtCore import QUrl 
from PyQt4.QtGui import QDesktopServices
import subprocess, os, guts

class MapMaker(Snuffling):
    '''
    <h1>Create a map containing event and station locations using googlemaps.</h1>

    Invokes the standard browser. Some browsers do not allow javascript to open
    and read the xml-file containing the necessary information due to the "Some-Origin-Policy".
    In that case you need to reset your standard browser.
    I.e.: Firefox on Linux do:
        xdg-settings set default-web-browser firefox.desktop
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
                active_event, active_stations = self.get_active_event_and_stations()
                break
            except AttributeError:
                active_stations = []
                active_event = None
                print 'Presumably no active event set'
                break

        station_list=[]
        stations= viewer.stations
        for NS, stat in stations.items():
            if '%s.%s.'%(str(NS[0]),str(NS[1])) in map(lambda x: x.nsl_string(), active_stations):
                xml_station_marker = XMLStationMarker(nsl = '%s.%s'%(str(NS[0]), str(NS[1])),
                                            longitude = stat.lon,
                                            latitude = stat.lat,
                                            active = 'yes')
            else:
                xml_station_marker = XMLStationMarker(nsl='%s.%s'%(str(NS[0]), str(NS[1])),
                                            longitude=stat.lon,
                                            latitude=stat.lat,
                                            active='no')
            station_list.append(xml_station_marker)
        active_station_list = StationMarkerList(stations=station_list)

        ev_marker_list = []
        if active_event is not None:
            xml_active_event_marker = XMLEventMarker(eventname = active_event.name,
                latitude=active_event.lat,
                longitude=active_event.lon,
                origintime=util.time_to_str(active_event.time),
                magnitude=active_event.magnitude,
                active='yes')
            ev_marker_list.append(xml_active_event_marker)    

        for m in markers:
            if isinstance(m, gui_util.EventMarker):
                ev = m.get_event()
                try:
                    xmleventmarker = XMLEventMarker(eventname=ev.name,
                                            longitude=ev.lon, 
                                            latitude=ev.lat, 
                                            origintime=util.time_to_str(ev.time), 
                                            magnitude=ev.magnitude,
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

        f=open(os.getenv("HOME")+'/.snufflings/googleMapsSnuffling/dumpedPyrockoMarker.xml', 'w')
        f.write(event_station_list.dump_xml())
        f.close()

        QDesktopServices.openUrl(QUrl('file://'+ os.getenv("HOME") + "/.snufflings/googleMapsSnuffling/readXml.html"))

def __snufflings__():
    return [ MapMaker() ]
