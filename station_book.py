from pyrocko.gui.snuffling import Snuffling, Choice


class StationBook(Snuffling):
    def setup(self):
        self.set_name("Station Book")
        self.attributes = [
            'network', 'station', 'location', 'lat', 'lon', 'elevation',
            'depth']
        self.add_parameter(
            Choice('Sort by', 'sort_by', 'network', self.attributes))

    def call(self):
        self.cleanup()

        stations = self.get_stations()
        stations = sorted(stations, key=lambda x: getattr(x, self.sort_by))
        ax = self.pylab()
        ax.axis('tight')
        ax.axis('off')
        cells = []

        col_labels = ('index', ) + tuple(self.attributes)
        row_format = "{:>10}" * (len(self.attributes) + 1)

        print(row_format.format(*col_labels))
        for i, s in enumerate(stations):
            d = [i]
            d.extend([getattr(s, attribute) for attribute in self.attributes])
            cells.append(d)
            print(row_format.format(*d))

        station_table = ax.table(
            cellText=cells,
            rowLabels=range(len(stations)),
            colLabels=col_labels,
            loc='center')

        for key, cell in station_table.get_celld().items():
                cell.set_linewidth(0)


def __snufflings__():
    return [StationBook()]
