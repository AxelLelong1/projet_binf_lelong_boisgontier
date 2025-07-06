import struct
import os
import sys
from asn1 import Decoder, Types

def read_null_terminated(data, offset):
    while data[offset] == 0:
        offset += 1
    start = offset
    while data[offset] != 0:
        offset += 1
    string = data[start:offset].decode('utf-8')
    return string, offset + 1

def read_sequence_type(nin_path):
    """
    Find the sequence and information
    """

    with open(nin_path, "rb") as f:
        data = f.read()

    # Read sequence type (4-byte integer at offset 4)
    format = struct.unpack_from(">I", data, 0)[0]

    sequence_type = struct.unpack_from(">I", data, 4)[0]

    offset = 12  # Skip format-version and sequence-type + 4 bits
    
    title, offset = read_null_terminated(data, offset)
    volume, offset = read_null_terminated(data, offset)
    create_date, offset = read_null_terminated(data, offset)
    while (offset % 4 != 0):
        offset += 1 # after create date, offset is aligned on base 4

    num_oids = struct.unpack_from('>I', data, offset)[0]
    offset += 4

    volume_length = struct.unpack_from('<Q', data, offset)[0]
    offset += 8

    max_seq_length = struct.unpack_from('>I', data, offset)[0]
    offset += 4

    print("Format Version:", format)
    print("Sequence type:", sequence_type)
    
    print("title:", title)
    print("volume name:", volume)
    print("create date: ", create_date)
    
    print("num_oids:", num_oids)
    print("volume_length:", volume_length)
    print("max_seq_length:", max_seq_length)
    
    print("Next byte (start of header-array):", data[offset])

    print("header 1:", struct.unpack_from('>I', data, offset)[0])
    print("header 2:", struct.unpack_from('>I', data, offset+4)[0])
    print("header 3:", struct.unpack_from('>I', data, offset+8)[0])

    offset+=num_oids

    print("Start of Sequence array:", data[offset])

    print("seq 1:", struct.unpack_from('>I', data, offset)[0])
    print("seq 2:", struct.unpack_from('>I', data, offset+4)[0])
    print("seq 3:", struct.unpack_from('>I', data, offset+8)[0])

    return num_oids


def main():
    if len(sys.argv) != 3:
        print("Usage: python dbread.py dbpath dbname")
        sys.exit(1)

    dbpath, dbname = sys.argv[1], sys.argv[2]

    nhr_path = os.path.join(dbpath, f"{dbname}.nhr")
    nin_path = os.path.join(dbpath, f"{dbname}.nin")
    nsq_path = os.path.join(dbpath, f"{dbname}.nsq")
    fasta_path = os.path.join(dbpath, f"{dbname}.fasta")

    num_oids = read_sequence_type(nin_path)


if __name__ == "__main__":
    main()