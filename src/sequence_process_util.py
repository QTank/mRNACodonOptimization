from hamiltonian import CodonOptimizer
from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon



import json


def get_qubit_len(seq, config, encode_type):
    if encode_type == "dense":
        return CodonOptimizer(seq, config, DenseCodon, "dense").qubit_len
    else:
        return CodonOptimizer(seq, config, OneHotCodon, "one-hot").qubit_len


def split_sequence_with_qubit_len(data, target_qubit_len, config, encode_type):
    chunks = []

    current = ""
    best_seq = ""
    best_diff = float("inf")

    for aa in data:
        current += aa

        q = get_qubit_len(current, config, encode_type)

        diff = abs(q - target_qubit_len)

        # 更新最优接近12的状态
        if diff < best_diff:
            best_diff = diff
            best_seq = current

        # 如果明显超出，就在最佳点切
        if q > target_qubit_len + 2:
            chunks.append(best_seq)
            print(best_seq, "qubit len:", CodonOptimizer(best_seq, config, DenseCodon, "dense").qubit_len)
            # 重置
            current = current[len(best_seq):]
            best_seq = ""
            best_diff = float("inf")

    # 收尾
    if current:
        chunks.append(current)

    return chunks


def test():
    with open("../config.json", "r") as f:
        config = json.load(f)

    data = "MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDV"
    res = split_sequence_with_qubit_len(data, 13, config, 'dense')
    print(len(res), res)

    for seq in res:
        print(seq, get_qubit_len(seq, config, 'one-hot'))
