from strenum import StrEnum
from enum import auto
from collections import deque
import paleobiodb_interface as pbdb
from paleobiodb_interface import int_vocab as rv

class TimeLevel(StrEnum):
    eon = auto()
    era = auto()
    period = auto()
    epoch = auto()
    # subepoch = 5
    age = auto()
    # subage = 7
    # zone = 8

    def index(self):
        cls = self.__class__
        members = list(cls)
        return members.index(self) + 1

    def next(self):
        '''Return the next lower level of time interval'''
        cls = self.__class__
        members = list(cls)
        index = members.index(self) + 1
        if index >= len(members):
            raise StopIteration('end of enumeration reached')
        return members[index]

    @classmethod
    def abbreviate_levels(cls, level):
        if level == cls.eon:
            return 'o'
        if level == cls.era:
            return 'r'
        if level == cls.period:
            return 'p'
        if level == cls.epoch:
            return 'e'
        if level == cls.age:
            return 's'

def queryColumn(search_lvl, progress_bar=None):
    # Initial query to get highest level intervals
    seedData = pbdb.query_geological_intervals()

    # Load intervals into stack LIFO (oldest on top)
    stack = deque()
    for record in seedData:
        if TimeLevel[record[rv.LEVEL]] == TimeLevel.eon:
            stack.append(record)

    def checkSubintervals(parent, childList):
        '''childList must be sorted youngest to oldest'''
        pointer = parent[rv.MIN_MA]
        for child in childList:
            if TimeLevel[child[rv.LEVEL]].index() != TimeLevel(parent[rv.LEVEL]).next().index() or child[rv.PARENT] != parent[rv.ID]:
                continue
            if child[rv.MIN_MA] != pointer:
                return False
            pointer = child[rv.MAX_MA]
        return parent[rv.MAX_MA] == pointer

    column = deque()
    while True:
        if len(stack) <= 0:
            break
        interval = stack.pop()

        if TimeLevel[interval[rv.LEVEL]].index() >= search_lvl.index():
            column.append(interval)
            if progress_bar is not None:
                progress_bar.update(1)
            continue

        subintervals = pbdb.query_geological_intervals(interval[rv.MIN_MA], interval[rv.MAX_MA])

        if checkSubintervals(interval, subintervals):
            for subint in subintervals:
                if TimeLevel[subint[rv.LEVEL]].index() == TimeLevel(interval[rv.LEVEL]).next().index():
                    stack.append(subint)
        else:
            column.append(interval)
            if progress_bar is not None:
                progress_bar.update(1)
    return column

def command_line_interface():
    import argparse
    import pickle
    from tqdm import tqdm

    expected_intervals = {TimeLevel.eon: 4,
                          TimeLevel.era: 11,
                          TimeLevel.period: 27,
                          TimeLevel.epoch: 49,
                          TimeLevel.age: 117}

    parser = argparse.ArgumentParser(description='Query PaleobioDB for geological column data.')
    parser.add_argument('-l', '--level', type=str, choices=[lvl.name for lvl in TimeLevel], default='age',
                        help='The maximum level of time intervals to include in the column (default: age)')
    parser.add_argument('-o', '--output', type=str, default='geological_column.pkl',
                        help='Output filename for the geological column data (default: geological_column.pkl)')
    parser.add_argument('--ls', action='store_true',
                        help='List time intervals and exit')
    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input filename for existing geological column data (default: None). If provided, the query will be skipped and data loaded from this file.')

    args = parser.parse_args()
    search_lvl = TimeLevel[args.level]

    if args.input is not None:
        with open(args.input, 'rb') as f:
            column = pickle.load(f)
        print(f'Geological column data loaded from {args.input}')
    else:
        print(f'Downloading geological column data from paleobiodb.org up to level: {search_lvl.name}...')
        t = tqdm(total=expected_intervals[search_lvl])
        column = queryColumn(search_lvl, progress_bar=t)
        t.close()

    if args.ls:
        print(f'Geological intervals up to level: {search_lvl.name}')
        for interval in column:
            print(f"{interval[rv.NAME]} ({interval[rv.MIN_MA]} Ma - {interval[rv.MAX_MA]} Ma) Level: {interval[rv.LEVEL]}")
    else:
        with open(args.output, 'wb') as f:
            pickle.dump(column, f)

        print(f'Geological column data saved to {args.output}')

if __name__ == '__main__':
    command_line_interface()