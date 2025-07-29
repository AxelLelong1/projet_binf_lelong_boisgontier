import struct
import os
import sys
from asn1 import Decoder, Types


# ---------- NIN FILE ---------- #

def read_len_prefixed_string(data, offset):
    """
    Lit une chaîne préfixée par sa longueur (1 octet), puis retourne la chaîne et l'offset suivant.
    """
    length = data[offset]
    offset += 1
    string = data[offset:offset + length].decode('utf-8')
    offset += length
    return string, offset

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
    while (data[offset] == 0):
        offset += 1
    title, offset = read_len_prefixed_string(data, offset)
    while (data[offset] == 0):
        offset += 1
    volume, offset = read_len_prefixed_string(data, offset)
    while (data[offset] == 0):
        offset += 1
    create_date, offset = read_len_prefixed_string(data, offset)
    while (offset % 4 != 0):
        offset += 1 # after create date, offset is aligned on base 4

    num_oids = struct.unpack_from('>I', data, offset)[0]
    offset += 4

    volume_length = struct.unpack_from('<Q', data, offset)[0]
    offset += 8

    max_seq_length = struct.unpack_from('>I', data, offset)[0]
    offset += 4

    """ print("Format Version:", format)
    print("Sequence type:", sequence_type)
    
    print("title:", title)
    print("volume name:", volume)
    print("create date: ", create_date)
    
    print("num_oids:", num_oids)
    print("volume_length:", volume_length)
    print("max_seq_length:", max_seq_length)
    
    print("Next byte (start of header-array):", offset)

    print("header 1:", struct.unpack_from('>I', data, offset)[0])
    print("header 2:", struct.unpack_from('>I', data, offset+4)[0])
    print("header 3:", struct.unpack_from('>I', data, offset+8)[0]) """

    offset+= (num_oids + 1) * 4

    # === Sequence Array ===

    """ print("Start of Sequence array:", offset)

    print("seq 1:", struct.unpack_from('>I', data, offset)[0])
    print("seq 2:", struct.unpack_from('>I', data, offset+4)[0])
    print("seq 3:", struct.unpack_from('>I', data, offset+8)[0]) """

    sequence_array = []
    for i in range(num_oids + 1):  # +1 pour le dernier offset
        value = struct.unpack_from('>I', data, offset)[0]
        sequence_array.append(value)
        offset += 4

    # === Ambig Array ===

    """ print("Start of Ambig Array:", offset)

    print("ambig 1:", struct.unpack_from('>I', data, offset)[0])
    print("ambig 2:", struct.unpack_from('>I', data, offset+4)[0])
    print("ambig 3:", struct.unpack_from('>I', data, offset+8)[0]) """

    ambig_array = []
    for i in range(num_oids + 1):  # idem
        value = struct.unpack_from('>I', data, offset)[0]
        ambig_array.append(value)
        offset += 4


    return sequence_array, ambig_array

# ---------- NSQ FILE ---------- #

def decode_ncbi_na2(byte_array, debug=False):
    """Decode NCBI Na2 packed nucleotide sequence (4 bases per byte).

    Args:
        byte_array: bytes object containing the packed sequence.

    Returns:
        String of nucleotides A,C,G,T.
    """
    if not byte_array:
        return ""

    # Last byte: last 2 bits indicate how many bases are valid in the last byte (0 means all 4 are valid)
    remainder = byte_array[-1] & 0b11
    packed_data = byte_array

    seq = ""
    for i, b in enumerate(packed_data):
        if b == 0 and i == len(packed_data) - 1:
            continue
        # For all bytes except the last, decode all 4 bases
        if i < len(packed_data) - 1:
            for shift in [6, 4, 2, 0]:
                bits = (b >> shift) & 0b11
                seq += "ACGT"[bits]
        else:
            # For last byte, only decode 'remainder' bases if remainder != 0, else decode all 4
            num_bases = remainder if remainder != 0 else 4
            for shift in [6, 4, 2, 0][:num_bases]:
                bits = (b >> shift) & 0b11
                seq += "ACGT"[bits]

    if debug:
        print(f"Decoded sequence: {seq}")
        print(f"Last byte remainder: {remainder}")
        print(f"Packed data length: {len(packed_data)} bytes")
        print(f"Packed data: {[f'{b:08b}' for b in packed_data]}")

    return seq

def read_sequences(seq, ambig, nsq_path):
    """
    Read the sequences
    """

    sequences = []

    with open(nsq_path, "rb") as f:
        
        for i in range(len(seq) - 1):
            start = int(seq[i])
            end = int(ambig[i])
            length_bytes = (end - start)

            f.seek(start)
            raw = f.read(length_bytes)
            sequences.append(decode_ncbi_na2(raw))
    
    return sequences

# ---------- NHR FILE ---------- # 

def read_ids_from_file(nhr_path):
    with open(nhr_path, "rb") as f:
        data = f.read()

    decoder = Decoder()
    decoder.start(data)
    return read_ids(decoder)

def find_accession_ids(obj):
    results = []
    if isinstance(obj, list):
        for item in obj:
            results.extend(find_accession_ids(item))
    elif isinstance(obj, bytes):
        try:
            s = obj.decode('ascii')
            if s and s[0].isalpha() and ' ' not in s:
                results.append(s)
        except Exception:
            pass
    return results

def read_ids(decoder, indent=0):
    ids = []
    while not decoder.eof():
        tag = decoder.peek()
        tag_class = tag.cls
        tag_nr = tag.nr

        #debug purpose: print('  ' * indent + f'Tag class: {tag_class}, number: {tag_nr}')

        if tag_class == Types.Constructed:
            decoder.enter()
            ids.extend(read_ids(decoder, indent + 1))  # recurse with same decoder
            decoder.leave()
        else:
            _, value = decoder.read()
            # debug purpose: print('  ' * indent + f'Value: {value}')
            ids.append(find_accession_ids(value))
    return ids

# ---------- FASTA FILE ---------- #

def print_fasta(sequences, ids, fasta_path):
    with open(fasta_path, "w") as f:
        for seq_id, seq in zip(ids, sequences):
            f.write(f">{seq_id[0]}\n")
            f.write(f"{seq}\n")
    return None

# ---------- MAIN FILE ---------- #

def main():
    if len(sys.argv) != 3:
        print("Usage: python dbread.py dbpath dbname")
        sys.exit(1)

    dbpath, dbname = sys.argv[1], sys.argv[2]

    nin_path = os.path.join(dbpath, f"{dbname}.nin")
    nsq_path = os.path.join(dbpath, f"{dbname}.nsq")
    nhr_path = os.path.join(dbpath, f"{dbname}.nhr")
    fasta_path = os.path.join(dbpath, f"{dbname}.fasta")

    seq, ambig = read_sequence_type(nin_path)
    sequences = read_sequences(seq, ambig, nsq_path)
    ids = read_ids_from_file(nhr_path)

    print_fasta(sequences, ids, fasta_path)

if __name__ == "__main__":
    main()