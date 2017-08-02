from pyrocko import model
from pyrocko.gui.snuffling import Snuffling


class ExtractEvents(Snuffling):
    '''
    Extract Events of selected Event Markers and write them to a catalog file.
    '''

    def setup(self):
        self.set_name('Write Events to Catalog')
        self.set_live_update(False)

    def call(self):
        markers = self.get_selected_event_markers()
        events = [m.get_event() for m in markers]
        if len(events) == 0:
            self.fail('no events found')

        out_filename = self.output_filename('Template for output files')
        model.dump_events(events, filename=out_filename)


def __snufflings__():
    return [ExtractEvents()]
