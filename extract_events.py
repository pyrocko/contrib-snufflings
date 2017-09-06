import sys
from pyrocko import util, io
from pyrocko.snuffling import Snuffling, load_markers, Param, NoViewerSet


default_output_filename = '%(eventname)s_%(network)s.%(station)s.%(location)s.%(channel)s.mseed'

class ExtractEvents(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h2 align="center">Extract Waveforms of Events</h2>
    <p>
    Select all event markers for which you want to extract waveforms. Set the
    <b>Start time rel. to event</b> and <b>End time rel. to event</b> and hit
    <b>Run</b>. You can modify the filename template to store
    chopped trace segments in the opening file dialog. Click <b>Save</b>.
    </p>
    </body>
    </html>
    '''

    def setup(self):
        self.set_name('Extract Events')
        self.add_parameter(Param('Start time rel. to event [s]', 'tbeg', 0., -3600., 0.))
        self.add_parameter(Param('End time rel. to event [s]', 'tend', 1200., 0., 3*3600.))
        self.set_live_update(False)

    def call(self):
        p = self.get_pile()

        try:
            markers = self.get_selected_event_markers()
        except NoViewerSet:
            markers = load_markers(self.markers_filename)

        try:

            out_filename = self.output_filename('Template for output files',
                                                default_output_filename)
        except NoViewerSet:
            out_filename = self.out_filename

        for m in markers:
            event = m.get_event()
            eventname = event.name
            if not eventname:
                eventname = util.time_to_str(event.time, format='%Y-%m-%d_%H-%M-%S')

            traces = p.all(tmin=event.time + self.tbeg,
                           tmax=event.time + self.tend)

            io.save(traces, out_filename, additional=dict(
                eventname=eventname))

    def configure_cli_parser(self, parser):

         parser.add_option(
            '--markers',
            dest='markers_filename',
            metavar='FILENAME',
            help='Read markers from FILENAME')

         parser.add_option(
            '--output',
            dest='out_filename',
            default=default_output_filename,
            metavar='FILENAME',
            help='set output filename template (default="%s")' % default_output_filename)


def __snufflings__():
    return [ExtractEvents()]


if __name__ == '__main__':
    import logging
    logger = logging.getLogger()
    util.setup_logging('extract_events.py', 'info')
    s = ExtractEvents()
    options, args, parser = s.setup_cli()
    s.markers_filename = options.markers_filename
    s.out_filename = options.out_filename
    if not options.markers_filename:
        logger.critical('no markers file given; use the --markers=FILENAME option')
        sys.exit(1)

    s.call()

