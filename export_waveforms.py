from pyrocko import io
from pyrocko.snuffling import Snuffling, NoViewerSet, Choice


class ExportWaveforms(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>

    <h1 align="center">Export selected or visible traces</h1>
    <body>
    <p>
    Choose the desired format from the <b>Format</b> menu and press
    <b>Run</b>.
    If no traces have been selected using extended markers all traces visible
    in the viewer will be exported.
    </p>
    <p>
    Note that exporting to miniseed requires the network, station, location and
    channel codes to be of length 2, 5, 2 and 3, respectively. Codes exceeding 
    these lenghts will be silently truncated.\br
    In order to have more control on code replacements it is recommended to use
    the command line tool <b>jackseis<b> which is shipped with pyrocko.
    </p>
    </body>
    </html>
    '''

    def setup(self):
        self.set_name('Export Waveforms')
        self.add_parameter(Choice('Format', 'format', 'mseed',
                                  ['mseed', 'text', 'sac', 'yaff']))
        self.set_live_update(False)

    def call(self):
        template = 'trace_%(network)s.%(station)s.%(location)s.%(channel)s'
        try:

            if self.format == 'text':
                default_output_filename = template + '.dat'

            else:
                default_output_filename = template + '.' + self.format

            out_filename = self.output_filename('Template for output files',
                                                default_output_filename)
        except NoViewerSet:
            out_filename = self.out_filename

        traces = self.chopper_selected_traces(fallback=True)
        for trs in traces:
            for tr in trs:
                if self.format == 'mseed':
                    if len(tr.network) > 2:
                        tr.set_network(tr.network[:2])
                    if len(tr.station) > 5:
                        tr.set_station(tr.station[:5])
                    if len(tr.location) > 2:
                        tr.set_location(tr.location[:2])
                    if len(tr.channel) > 3:
                        tr.set_channel(tr.channel[:3])
                io.save(tr, out_filename, format=self.format)

def __snufflings__():
    return [ExportWaveforms()]


