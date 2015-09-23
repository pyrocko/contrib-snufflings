from pyrocko.snuffling import Snuffling, Param
import sys 
import vtk 
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as num 
from pyrocko import model
from pyrocko import orthodrome as ortho
import os
from vtk_focmec import *

class VtkTest(Snuffling):
    def setup(self):
        self.set_name('VTK')
        self.add_parameter(Param('Map Radius [km]', 'radius', 10, 10, 100))

    def call(self):
        self.cleanup()
        viewer = self.get_viewer()
        markers = self.get_selected_event_markers()
        events = [m.get_event() for m in markers]
        lats = [e.lat for e in events]
        lons = [e.lon for e in events]
        center_lat = num.mean(lats)
        center_lon = num.mean(lons)
        #tmin, tmax = self.get_selected_time_range(fallback=True)
        #event_markers = filter(lambda x: x.tmin>=tmin and x.tmax<=tmax,
        #                       viewer.markers)

        #event_markers = filter(lambda x: isinstance(x, EventMarker), event_markers)
        #if self.maxd:
        #    event_markers = filter(lambda x: distance(self, x._event)<=self.maxd*km,
        #                       event_markers)

        #if event_markers==[]:
        #    self.fail('No events in selected area found')

        self.frame = QtGui.QFrame(self._panel_parent)

        self.vl = QtGui.QGridLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vtkWidget.setMinimumSize(900, 880)
        self.vl.addWidget(self.vtkWidget)

        self.ren = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        actors = []
        normals_list, centers, colors = read_data(events=events, get=['p-axis'])
        for i, normals in enumerate(normals_list):
            kwargs = {"centers": centers, 'normals':normals, "return_pdm":True}

        colors = [(0,1,0)]*len(colors)
        kwargs = {"centers": centers, 'normals':normals, "return_pdm":True}
        actor1, apd = make_polydata_actor(**kwargs)
        actors.append(actor1)
        map_actor = setup_vtk_map_actor(center_lat, center_lon, self.radius * 1000.)
        actors.append(map_actor)
        setup_renderer(self.ren, actors, bboxpolydata=apd)

        self.ren.ResetCamera()

        #self.ren.SetViewport(0.,0.,1.,1.)
        self._panel_parent.add_tab('pressure', self.frame)

        self.iren.Initialize()

def __snufflings__():
    return [VtkTest()]

