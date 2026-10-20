#!/usr/bin/env python3
"""Write a GollumFit compact (FastMC) file whose simulated events are the input's repeated k times.

Used only for memory-scaling tests: physics is unchanged up to an overall factor k in the
expectation. Format (src/compactIO.cpp): uint32 program checksum, uint64 n_data, n_data x Event,
uint64 n_sim, n_sim x Event, uint32 CRC-32 of everything before it. sizeof(Event) = 344 here.
Usage: python3 make_scaled_compact.py <in.fastmc> <out.fastmc> <k>
"""
import os, struct, sys, zlib

src, dst, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
size = os.path.getsize(src)
with open(src, "rb") as f:
    prog, nexp = struct.unpack("<IQ", f.read(12))
    data_bytes = f.read(nexp * 344)
    nsim = struct.unpack("<Q", f.read(8))[0]
    ev_start = f.tell()
    assert ev_start + nsim * 344 + 4 == size, "unexpected layout"
with open(src, "rb") as f, open(dst, "wb") as out:
    crc = 0
    def w(b):
        global crc
        out.write(b); crc = zlib.crc32(b, crc)
    w(struct.pack("<IQ", prog, nexp)); w(data_bytes); w(struct.pack("<Q", nsim * k))
    for _ in range(k):
        f.seek(ev_start); remaining = nsim * 344
        while remaining:
            chunk = f.read(min(1 << 24, remaining)); w(chunk); remaining -= len(chunk)
    out.write(struct.pack("<I", crc))
print("wrote %s: %d simulated events (%.2f GB)" % (dst, nsim * k, os.path.getsize(dst) / 1e9))
