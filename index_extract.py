"""
该脚本用于提取illumina nextseq500/550测序数据中的接头序列，最终生成一个index.txt文件，包含接头以及其出现次数

"""
import gzip
import fire
import struct
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


class BCI:
    def __init__(self, file):
        self.file = Path(file)

    def parse(self):
        tiles = []
        with open(self.file, 'rb') as f:
            data = f.read()
        offset = 0
        tile_count = len(data) // 8
        print(f'tile count: {tile_count}')
        for x in range(tile_count):
            tid, count = struct.unpack('II', data[x*8:x*8+8])
            tiles.append(Tile(tid, count, offset))
            offset += count
        return tiles


class Tile:
    def __init__(self, tid, count, offset):
        self.tid = tid
        self.count = count
        self.offset = offset
        self.cycle_datas = []

    def load_bcl(self, bcls):
        for bcl in bcls:
            self.cycle_datas.append(bcl.get_tile_data(self))

    def extract(self):
        print(f"[{self.tid}] extract...")
        m = ['A', 'C', "G", 'T']
        seqs = []
        for i in range(len(self.cycle_datas[0])):
            seq = []
            for buff in self.cycle_datas:
                b = buff[i]
                if b == 0:
                    seqs.append('N')
                    continue
                base = m[b & 0b00000011]
                seq.append(base)
            seqs.append(''.join(seq))
        return seqs


class BCL:
    def __init__(self, file, bci: BCI):
        self.file = Path(file)
        print(f"loading {self.file}...")
        with gzip.open(self.file, 'rb') as f:
            self.data = f.read()
        self.header = self.data[:4]
        self.body = self.data[4:]
        self.bci = bci

    def get_tile_data(self, tile):
        return bytearray(self.body[tile.offset:tile.offset+tile.count])


def extract_lane(lane_dir, start, length, outfile):
    lane_num = lane_dir[-1]
    lane_dir = Path(lane_dir)
    bci = BCI(lane_dir / f's_{lane_num}.bci')
    tiles = bci.parse()
    print(f"tile count: {len(tiles)}")
    bcls = []
    for cycle in range(start, start+length):
        cycle_file = lane_dir / f'{cycle:04}.bcl.bgzf'
        bcl = BCL(cycle_file, bci)
        bcls.append(bcl)
    for tile in tiles:
        tile.load_bcl(bcls)

    print('start parse')
    results = []
    with ProcessPoolExecutor(max_workers=16) as ex:
        for tile in tiles:
            results.append(ex.submit(tile.extract))

    counter = Counter()
    for future in results:
        res = future.result()
        counter.load(res)
    counter.write(outfile)


class Counter:
    def __init__(self):
        self.counter = defaultdict(int)

    def load(self, seqs):
        for seq in seqs:
            self.counter[seq] += 1

    def write(self, outfile):
        with open(outfile, 'w') as f:
            for seq, count in self.counter.items():
                f.write(f'{seq[:8]}+{seq[8:]}\t{count}\n')


if __name__ == '__main__':
    fire.Fire(extract_lane)
