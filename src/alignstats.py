import sys
import math

def parse_args():
    if len(sys.argv) != 7:
        print("Usage: python alignstats.py Nseq Nnuc seq rwd,pen score cmd")
        sys.exit(1)
    Nseq = int(sys.argv[1])
    Nnuc = int(sys.argv[2])
    seq = sys.argv[3].upper()
    rwd, pen = map(int, sys.argv[4].split(','))
    score = int(sys.argv[5])
    cmd = sys.argv[6]
    return Nseq, Nnuc, seq, rwd, pen, score, cmd

def count_freqs(seq):
    total = len(seq)
    freqs = {base: seq.count(base)/total for base in "ATGC"}
    return freqs

def get_score_matrix(rwd, pen):
    return {
        (b1, b2): (rwd if b1 == b2 else pen)
        for b1 in "ATGC" for b2 in "ATGC"
    }

def get_score_dist(freqs_seq, rwd, pen):
    fdb = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}
    sigma = get_score_matrix(rwd, pen)
    C = {}
    for b1 in "ATGC":
        for b2 in "ATGC":
            score = sigma[(b1, b2)]
            C[score] = C.get(score, 0) + freqs_seq[b1] * fdb[b2]
    # Ajouter les scores manquants
    for s in [1, 0, -1, -2]:
        if s not in C:
            C[s] = 0.0
    total = sum(C.values())
    if total == 0:
        P = {k: 0.0 for k in C}
    else:
        P = {s: c / total for s, c in C.items()}
    score_probs = [(s, P[s]) for s in sorted(P.keys(), reverse=True)]
    return score_probs, P

def f_poly(x, score_probs):
    return sum((p if s != 0 else p - 1) * (x ** i) for i, (s, p) in enumerate(score_probs))

def df_poly(x, score_probs):
    return sum(i * (p if s != 0 else p - 1) * (x ** (i - 1)) for i, (s, p) in enumerate(score_probs) if i > 0)

def solve_lambda(score_probs):
    tolx = 1e-5
    x = math.exp(-0.5)
    a, b = 0.0, 1.0
    isNewt = 0
    f = 99
    for k in range(20):
        fprev = f
        f = f_poly(x, score_probs)
        df = df_poly(x, score_probs)

        if f > 0:
            a = x
        elif f < 0:
            b = x
        else:
            break

        if b - a < 2 * a * (1 - b) * tolx:
            x = (a + b) / 2
            break

        if (isNewt and abs(f) >= 0.9 * abs(fprev)) or df >= 0:
            x = (a + b) / 2
        else:
            p = -f / df
            y = x + p
            if y <= a or y >= b:
                x = (a + b) / 2
            else:
                isNewt = 1
                x = y
                if abs(p) < tolx * x * (1 - x):
                    break
    return -math.log(x)

def compute_H(lmbda, smax, score_probs):
    total = sum(s * p * math.exp(-i * lmbda) for i, (s, p) in enumerate(score_probs))
    denom = math.exp(-lmbda * smax)
    return (lmbda / denom) * total if denom != 0 else 0.0

def compute_K(H, lmbda, sigma, smax):
    if H == 0 or lmbda == 0:
        return 0.0
    if smax == 1:
        return H / lmbda * (1 - math.exp(-lmbda))
    avg_sigma = sum(sigma[(b1, b2)] for b1 in "ATGC" for b2 in "ATGC") / 16
    return (avg_sigma ** 2) * lmbda / H * (1 - math.exp(-lmbda))

def compute_alpha(H, lmbda):
    return lmbda / H if H != 0 else 0.0

def compute_beta(smin, smax):
    # Cas classiques d'alignement local où beta vaut -2
    if (smax == 1 and smin == -1) or (smax == 2 and smin == -3):
        return -2
    # Cas où smax=1 et smin=-2 (comme souvent 1,-2)
    if smax == 1 and smin == -2:
        return -2
    # Autres cas à vérifier si besoin, sinon 0
    return 0

def compute_adj(Nseq, Nnuc, M, alpha, beta, lmbda, K):
    mb = M * Nseq + Nnuc
    c = Nnuc * M - max(Nnuc, M) / K

    # Discriminant pour racine carrée (protection)
    discriminant = mb ** 2 - 4 * Nseq * c
    if discriminant < 0:
        return 1  # valeur de sécurité

    emax = 2 * c / (mb + math.sqrt(discriminant))
    emin = 0

    def h(e):
        val = (M - e) * (Nnuc - Nseq * e)
        if val <= 0:
            # log impossible => retourner un grand négatif pour forcer e < h(e)
            return float('-inf')
        return beta + (alpha / lmbda) * (math.log(K) + math.log(val))

    # dichotomie classique
    for _ in range(30):
        mid = (emin + emax) / 2
        hmid = h(mid)
        if hmid > mid:
            emin = mid
        else:
            emax = mid
        if abs(hmid - mid) < 1e-5:
            break

    res = int(math.floor(emin))

    # Vérification complémentaire (optionnelle)
    e_ceil = math.ceil(emin)
    if e_ceil <= emax and h(e_ceil) >= e_ceil:
        res = int(math.floor(e_ceil))

    return res

def compute_searchsp(Nseq, Nnuc, M, alpha, beta, lmbda, K):
    adj = compute_adj(Nseq, Nnuc, M, alpha, beta, lmbda, K)
    return int((Nnuc - Nseq * adj) * (M - adj))

def compute_E(Nsearch, lmbda, K, score):
    if K <= 0:
        return float('inf')
    return Nsearch * math.exp(-lmbda * score + math.log(K))

def compute_bitscore(score, lmbda, K):
    if K <= 0:
        return 0
    return int((score * lmbda - math.log(K)) / math.log(2))

def main():
    Nseq, Nnuc, seq, rwd, pen, score, cmd = parse_args()
    freqs = count_freqs(seq)
    score_probs_list, score_probs_dict = get_score_dist(freqs, rwd, pen)
    smax = max(score_probs_dict.keys())
    smin = min(score_probs_dict.keys())
    lmbda = solve_lambda(score_probs_list)
    H = compute_H(lmbda, smax, score_probs_list)
    sigma = get_score_matrix(rwd, pen)
    K = compute_K(H, lmbda, sigma, smax)
    alpha = compute_alpha(H, lmbda)
    beta = compute_beta(smin, smax)
    M = len(seq)
    searchsp = compute_searchsp(Nseq, Nnuc, M, alpha, beta, lmbda, K)
    Evalue = compute_E(searchsp, lmbda, K, score)
    bitscore = compute_bitscore(score, lmbda, K)

    if cmd == "lambda":
        print(f"{lmbda:.2f}")
    elif cmd == "H":
        print(f"{H:.2f}")
    elif cmd == "K":
        print(f"{K:.2f}")
    elif cmd == "alpha":
        print(f"{alpha:.2f}")
    elif cmd == "beta":
        print(f"{beta:.2f}")
    elif cmd == "searchsp":
        print(f"{searchsp}")
    elif cmd == "E":
        print(f"{Evalue:.0e}")
    elif cmd == "bitscore":
        print(f"{bitscore}")
    elif cmd == "all":
        print(f"{lmbda:.2f} {H:.2f} {K:.2f} {alpha:.2f} {beta:.2f} {searchsp} {Evalue:.0e} {bitscore}")
    else:
        print("Invalid cmd")

if __name__ == "__main__":
    main()
