import os
import sys
import tempfile
from dbread import read_sequence_type, read_sequences, read_ids_from_file, print_fasta
from dbalign import generate_hsw, fasta_to_dict, generate_hsp, extend_hsp
from alignstats import get_score_dist, solve_lambda, compute_H, get_score_matrix, compute_K, compute_alpha, compute_beta, compute_searchsp, compute_E, compute_bitscore, count_freqs

def main():
    if len(sys.argv) != 4:
        print("Usage: python myblast.py dbpath dbname query_seq")
        sys.exit(1)

    dbpath, dbname, query_seq = sys.argv[1], sys.argv[2], sys.argv[3].upper()
    nin = os.path.join(dbpath, dbname + ".nin")
    nsq = os.path.join(dbpath, dbname + ".nsq")
    nhr = os.path.join(dbpath, dbname + ".nhr")

    # --- Étape 1 : Lire la base NCBI ---
    seq_offsets, ambig = read_sequence_type(nin)
    sequences = read_sequences(seq_offsets, ambig, nsq)
    ids = read_ids_from_file(nhr)

    # --- Étape 2 : Écriture temporaire en FASTA ---
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmpf:
        fasta_path = tmpf.name
        print_fasta(sequences, ids, fasta_path)

    # --- Étape 3 : Alignement ---
    W = 11
    hsw = generate_hsw(query_seq, W)
    fasta_dict = fasta_to_dict(fasta_path)
    hsp = generate_hsp(hsw, fasta_dict, W)

    results = []

    for subject_id, subject_pos, query_pos, word in hsp:
        res = extend_hsp(query_seq, fasta_dict[subject_id], query_pos, subject_pos, W)
        if res is not None:
            q_start, q_end, s_start, s_end, taille, score, align_q, align_s = res

            # --- Étape 4 : Calcul E-value ---
            Nseq = len(fasta_dict)
            Nnuc = sum(len(seq) for seq in fasta_dict.values())
            M = len(query_seq)
            rwd, pen = 1, -2

            freqs = count_freqs(query_seq)
            score_probs, prob_dict = get_score_dist(freqs, rwd, pen)
            smax = max(prob_dict.keys())
            smin = min(prob_dict.keys())
            sigma = get_score_matrix(rwd, pen)
            lmbda = solve_lambda(score_probs)
            H = compute_H(lmbda, smax, score_probs)
            K = compute_K(H, lmbda, sigma, smax)
            alpha = compute_alpha(H, lmbda)
            beta = compute_beta(smin, smax)
            searchsp = compute_searchsp(Nseq, Nnuc, M, alpha, beta, lmbda, K)
            Evalue = compute_E(searchsp, lmbda, K, score)
            bitscore = compute_bitscore(score, lmbda, K)

            results.append({
                "subject_id": subject_id,
                "s_start": s_start,
                "s_end": s_end,
                "q_start": q_start,
                "q_end": q_end,
                "taille": taille,
                "score": score,
                "Evalue": Evalue,
                "bitscore": bitscore,
                "lambda": lmbda,
                "H": H,
                "K": K,
                "alpha": alpha,
                "beta": beta,
                "searchsp": searchsp,
                "align_q": align_q,
                "align_s": align_s
            })

    # --- Étape 5 : Filtrage par E-value minimale ---
    if not results:
        print("Aucun alignement significatif trouvé.")
        return

    min_eval = min(r["Evalue"] for r in results)
    kept = [r for r in results if abs(r["Evalue"] - min_eval) < 1e-15]

    for r in kept:
        print(f"{r['subject_id']} {r['s_start']} {r['s_end']} {r['q_start']} {r['q_end']} {r['taille']} {r['score']}")
        print(f"{r['lambda']:.2f} {r['H']:.2f} {r['K']:.2f} {r['alpha']:.2f} {r['beta']:.2f} {r['searchsp']} {r['Evalue']:.0e} {r['bitscore']}")
        print(r['align_q'])
        print(r['align_s'])

    # Nettoyage du fichier temporaire
    os.remove(fasta_path)

if __name__ == "__main__":
    main()
