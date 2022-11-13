from pyrocko.gui.snuffling import Snuffling, Choice, Param
from pyrocko.gui.marker import PhaseMarker
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
			None, 'positive', 'negative', 'undecidable']))

		self.set_live_update(True)

	def call(self):
		'''Main work routine of the snuffling.'''

		viewer = self.get_viewer()

		if viewer.selected_markers() == []:
			self.error('Error: No marker selected!')
		elif isinstance(viewer.selected_markers()[0], PhaseMarker) == False:
			self.error('Error: Selected marker is no phase marker!')

		markers = [
			m for m in viewer.selected_markers()
			if isinstance(m, PhaseMarker)]
		for m in markers:
			m.set_phasename(self.phase)
			if self.polarity == 'positive':
				m.set_polarity(1)
			elif self.polarity == 'negative':
				m.set_polarity(-1)
			elif self.polarity == 'undecidable':
				m.set_polarity(0)
			else:
				m.set_polarity(self.polarity)
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
