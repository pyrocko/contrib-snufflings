from pyrocko.gui.snuffling import Snuffling, Choice, Param
import pyrocko.gui.marker
#from pyrocko.pile_viewer import Marker


class Picker(Snuffling):

	'''
	Snuffling for having a panel to control the picking procedure, in
	addition to the key-based control.
	'''

	def setup(self):
		'''Customization of the snuffling.'''

		self.set_name('Picker')
		self.add_parameter(Choice('Phase', 'phase', None, [
			None, 'P', 'S', 'Pn', 'Pg', 'Sn', 'Sg']))
		self.add_parameter(Choice('Category', 'category', 'automatic', [
			'automatic', 'manual', 'revised', 'abandoned']))
		self.add_parameter(Choice('Polarity', 'polarity', None, [
			None, '+', '-']))
		self.add_parameter(Param(
			'pick adjustment [s]', 'tcorrect', 0., -5., 5.))
		self.first = True

		self.set_live_update(True)

	def call(self):
		'''Main work routine of the snuffling.'''

		self.cleanup()
		if self.phase is not None:
			Phase = self.phase

		viewer = self.get_viewer()

		# if viewer.get_selected_markers == [] and isinstance(viewer.get_selected_markers(), pyrocko.gui.marker.PhaseMarker):
		# 	print('Error: No marker selected!')

		markers = [
			m for m in viewer.selected_markers()
			if isinstance(m, pyrocko.gui.marker.PhaseMarker)]
		for m in markers:
			if self.first:
				time = m.get_tmin()
				self.first = False
			m.set_phasename(self.phase)
			m.set_polarity(self.polarity)
			attr = m.get_attributes()
			print(m.get_tmin())
			# we have the problem, that tmin and tmax will be updated continously
			# either the viewer will be updated as well or tmin and tmax are static
			m.set(
				m.get_nslc_ids(),
				time + self.tcorrect,
				time + self.tcorrect)

			if self.category == 'automatic':
				m.set_kind(0)
			elif self.category == 'manual':
				m.set_kind(1)
			elif self.category == 'revised':
				m.set_kind(2)
			elif self.category == 'abandoned':
				m.set_kind(3)


def __snufflings__():
	'''Returns a list of snufflings to be exported by this module.'''

	return [Picker()]
