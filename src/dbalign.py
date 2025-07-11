import sys

# Matrice de substitution pour A,T,G,C
sub_matrix = {
    ('A', 'A'): 1, ('A', 'T'): -2, ('A', 'G'): -2, ('A', 'C'): -2,
    ('T', 'A'): -2, ('T', 'T'): 1, ('T', 'G'): -2, ('T', 'C'): -2,
    ('G', 'A'): -2, ('G', 'T'): -2, ('G', 'G'): 1, ('G', 'C'): -2,
    ('C', 'A'): -2, ('C', 'T'): -2, ('C', 'G'): -2, ('C', 'C'): 1,
}

# ---------- FUNCTIONS ---------- #

def generate_hsw(seq, W):
    hsw = dict()

    for i in range(len(seq) - W + 1):
        word = seq[i:i+W]
        if word not in hsw:
            hsw[word] = []
        hsw[word].append(i)
    return hsw

def fasta_to_dict(fasta_path):
    sequences = dict()

    with open(fasta_path) as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            s = lines[i].strip()
            if s[0] == '>':
                key = s[1:]
            else:
                sequences[key] = s
    return sequences

def generate_hsp(hsw, fasta_dict, W):
    hsp = []
    for id in fasta_dict:
        sequence = fasta_dict[id]
        for i in range(len(sequence) - W + 1):
            word = sequence[i:i+W]
            if word in hsw:
                for qpos in hsw[word]:
                    hsp.append((id, i, qpos, word))
    return hsp

def extend_hsp(seq_query, seq_subject, q_start, s_start, W):
    T = 5
    dropoff = 20
    score = sum(sub_matrix[(seq_query[q_start + k], seq_subject[s_start + k])] for k in range(W))
    best_score = score

    # Étendre à gauche
    i, j = q_start - 1, s_start - 1
    ql, sl = q_start, s_start
    current_score = score

    while i >= 0 and j >= 0:
        current_score += sub_matrix[(seq_query[i], seq_subject[j])]
        if current_score > best_score:
            best_score = current_score
            ql, sl = i, j
        elif current_score < best_score - dropoff or current_score <= T:
            break
        i -= 1
        j -= 1

    # Étendre à droite
    i, j = q_start + W, s_start + W
    qr, sr = q_start + W, s_start + W
    current_score = score

    while i < len(seq_query) and j < len(seq_subject):
        current_score += sub_matrix[(seq_query[i], seq_subject[j])]
        if current_score > best_score:
            best_score = current_score
            qr, sr = i + 1, j + 1
        elif current_score < best_score - dropoff or current_score <= T:
            break
        i += 1
        j += 1

    # Calcul final du score
    final_score = 0
    for k in range(qr - ql):
        final_score += sub_matrix[(seq_query[ql + k], seq_subject[sl + k])]

    if final_score > T:
        return (ql, qr-1, sl, sr-1, qr - ql, final_score,
                seq_query[ql:qr], seq_subject[sl:sr])
    return None


# ---------- MAIN PART ---------- #

def main():
    if len(sys.argv) != 3:
        print("Usage: python dbalign.py dbfasta seq")
        sys.exit(1)

    dbfasta = sys.argv[1]
    query_seq = sys.argv[2]
    W = 11  # Longueur des HSW

    hsw = generate_hsw(query_seq, W)
    fasta_dict = fasta_to_dict(dbfasta)
    hsp = generate_hsp(hsw, fasta_dict, W)

    seen = set()

    for subject_id, subject_pos, query_pos, word in hsp:
        result = extend_hsp(query_seq, fasta_dict[subject_id], query_pos, subject_pos, W)
        if result is not None:
            q_start, q_end, s_start, s_end, taille, score, align_q, align_s = result
            align_key = (subject_id, s_start, s_end, q_start, q_end, score)
            if align_key not in seen:
                seen.add(align_key)
                print(f"{subject_id} {s_start} {s_end} {q_start} {q_end} {taille} {score}")
                print(align_q)
                print(align_s)

# ------------ TESTS ------------ #
#    # HSW
#    hsw = generate_hsw(query_seq, W)
#    print("HSW:")
#    for word, positions in hsw.items():
#        print(f"{word}: {positions}")
#    print()
#
#    # HSP
#    fasta_dict = fasta_to_dict(dbfasta)
#    hsp = generate_hsp(hsw, fasta_dict, W)
#    print("HSP:")
#    for subject_id, subject_pos, query_pos, word in hsp:
#        print(f"{word} | DB: {subject_id}:{subject_pos} <-> HSW:{query_pos}")

if __name__ == "__main__":
    main()