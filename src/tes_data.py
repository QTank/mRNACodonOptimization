import util
from hamiltonian import CodonOptimizer
from denseCodon import DenseCodon

from oneHotCodon import OneHotCodon
import json

with open("../config.json", "r") as f:
    config = json.load(f)

time = [20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
time_qaoa = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4098, 8196]


def test(file_name, encode_type, chunk=7):
    data = util.split_sequence_from_file(file_name, chunk)
    print(len(util.parse_sequence_from_file(file_name)), data)
    min_len, max_len = 100, 0

    time_est, i = 0, 0
    for sequence in data:
        i += 1
        if encode_type == 'dense':
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        else:
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_len = codon_opt.qubit_len
        min_len = min(min_len, qubit_len)
        max_len = max(max_len, qubit_len)
        print(f"{i}/{len(data)} {sequence}: {qubit_len}, {CodonOptimizer(sequence, config, DenseCodon, 'dense').qubit_len}")
        if chunk <= 6:
            time_est += time_qaoa[qubit_len-5]
    print(f"min: {min_len}, max: {max_len}")
    print(f"estimate time: {time_est}")

def test_spilt():
    target = 12
    data = "SWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQV"
    data = "MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDV"
    chunks = []
    current_seq = ""
    prev_seq = ""
    prev_qubit = 0
    time_execute = 0
    for aa in data:
        current_seq += aa

        codon_opt = CodonOptimizer(current_seq, config, DenseCodon, "dense")
        qubit_len = codon_opt.qubit_len

        if qubit_len >= target:
            # choose closest to target
            if abs(prev_qubit - target) < abs(qubit_len - target):
                chunks.append(prev_seq)
                print(f"{prev_seq}, qubit len: {prev_qubit}, time: {time[prev_qubit - 11]}s")
                current_seq = aa
                time_execute += time[prev_qubit - 11]
            else:
                chunks.append(current_seq)
                print(f"{current_seq}, qubit len: {qubit_len}, time: {time[qubit_len - 11]}s")
                current_seq = ""
                time_execute += time[qubit_len - 11]

        prev_seq = current_seq
        prev_qubit = qubit_len

    print(f"time: {time_execute}")


def test_clean():
    target = 12
    chunks = []
    data = "MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDV"

    current = ""
    best_seq = ""
    best_diff = float("inf")
    time_execute = 0
    for aa in data:
        current += aa
        opt = CodonOptimizer(current, config, DenseCodon, "dense")
        q = opt.qubit_len

        diff = abs(q - target)

        # 更新最优接近12的状态
        if diff < best_diff:
            best_diff = diff
            best_seq = current

        # 如果明显超出，就在最佳点切
        if q > target + 2:
            chunks.append(best_seq)
            qubit_len = CodonOptimizer(best_seq, config, DenseCodon, "dense").qubit_len
            print(best_seq, "qubit len:", qubit_len)

            # 重置
            current = current[len(best_seq):]
            best_seq = ""
            best_diff = float("inf")
            time_execute += time[qubit_len - 11]
    # 收尾
    if current:
        chunks.append(current)

    print(f"time: {time_execute}")


file_name_log = ['01-sars2_spike_vaccine.fasta', '02-sars2_n_vaccine.fasta', '03-influenza_ha_vaccine.fasta',
                 '04-Zika_E_ectodomain.fasta', '05-DENV1_E_ectodomain.fasta', '06-rabies_vaccine_antigen.fasta',
                 '07-ebola_vaccine_antigen.fasta', '08-mers_vaccine_antigen.fasta', '09-nipah_vaccine_antigen.fasta',
                 '10-hendra_vaccine_antigen.fasta', '11-marburg_vaccine_antigen.fasta',
                 '12-lassa_vaccine_antigen.fasta',
                 '13-Poliovirus_VP1.fasta', '14-Norovirus_VP1_mRNA_vaccine.fasta', '15-HBV_HBsAg_vaccine.fasta']
for file_name in file_name_log:
    if file_name != '03-influenza_ha_vaccine.fasta': continue
    print(f'{file_name}')
    test("../data/" + file_name, 'one-hot', 3)
