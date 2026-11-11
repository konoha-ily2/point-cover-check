# Fetch the point-weight certificates by other authors, check their SHA256, and convert them into the format
# read by interval_cover.py. Also build the negative controls (certificates broken on purpose).
# Usage: python fetch_and_prep.py <work folder>
#   Original files (JSON converted to text) go to <work folder>/orig/, converted files to <work folder>/.
#   Then: python run_batch.py <work folder> --workers 2 --minutes 8
# The downloaded files are other authors' data and are not included in this repository.
import os, sys, hashlib, subprocess, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EV = "https://raw.githubusercontent.com/evand/square-packing/167d842cd27ba1451cb2833773ea930c80b9e65b/s12/certificates/"
SP = "https://raw.githubusercontent.com/squarepacker/s12-lower-bound/7a96bec36bc6811c3715ef581598f22ff9b7ba3a/"
WB = "https://raw.githubusercontent.com/wand125/square-packing-bounds/678aa1d029423177792f51d3e7e42895b749241e/certificates/"
FILES = [  # (saved name, URL, SHA256)
    ("s12_56points_3.8.txt", EV + "s12_56points_3.8.txt", "f78ac17e7a424d04952aa997b16ee5873e23fd1eb5e1bd1df641bde11d8abc83"),
    ("s12_uniform_7of81_3.888.txt", EV + "s12_uniform_7of81_3.888.txt", "728a681fb90a81bec8a2db2a482b7c47226453996de2023a937423d0f3000ce7"),
    ("s12_lower_3.931795_sparse.txt", EV + "s12_lower_3.931795_sparse.txt", "cc712a39643e124d63fe2be0f6763c105c5be1fdcbe9dfbfec8b9598e2726906"),
    ("s12_lower_3.931795.txt", EV + "s12_lower_3.931795.txt", "66ed8a0eba55e1be2d543982136b440d82af7fdd12ea0ed2fda521e4a4333d11"),
    ("s12_lower_3.9676.txt", EV + "s12_lower_3.9676.txt", "55c9422e6f01c0f3101ffba39135b2cc21c07a86e03e791e13afadd43cc13dd4"),
    ("s11_lower_3.8143.txt", EV + "s11_lower_3.8143.txt", "0b716e4079b07d333a42412021c970a97c861b933a76261a26cfe63841c55ac9"),
    ("s12_lower_3.9686.txt", EV + "s12_lower_3.9686.txt", "75f1cc891a8b8739b92b50c7ddbd85a493efa58cf2e921d58ce4ac7c2dcabe78"),
    ("s21_lower_4.9950.txt", EV + "s21/s21_lower_4.9950.txt", "c8e8f878205f2da9c213e4a87c06f17c5a759dce7e695e676dd5c50e0994f2ef"),
    ("s12_lower_3.969118.txt", SP + "s12_lower_3.969118.txt", "6ad9b0e8257166687f4e861b97f024f4993d7b2125897173b2f5e9f6e2167578"),
    ("s12_lower_3.9715.txt", SP + "s12_lower_3.9715.txt", "e2f326b28142cf22402f88357f4c7fe4680a08a32ae785b493adac335bcc2685"),
    ("56points_critical_760_199.txt", SP + "controls/56points_critical_760_199.txt", "d98681c28374e6af0febc7439dabef4bb0e4dae51921b5631002cdba63cfbaf4"),
    ("56points_over_1520_397.txt", SP + "controls/56points_over_1520_397.txt", "0b712361b78bcfe73a31ed409f643890596caeaca8cfedbb56e804e9e327f4fb"),
    ("scaled_further_31360_7900.txt", SP + "controls/scaled_further_31360_7900.txt", "04d7105c5b6d8b85cd0655f544b5c59163f9a36fa86814601a30c4114d1526f1"),
    ("control1_lowered_orbit.txt", SP + "controls/3.9715/control1_lowered_orbit.txt", "1c42bbb4cf21ee75d9e13e3947cde199559f5207f71662394ed66eecfd05aad9"),
    ("control2_scaled_further.txt", SP + "controls/3.9715/control2_scaled_further.txt", "f1796e271c6b877692ee06520153636cb7386e331724a720e966ec49c42b8560"),
    ("cert_n26_L545.json", WB + "cert_n26_L545.json", "9c7689036c3ec2c34e0038d81d9e2a08a77c0a2ff1c020991f1a5e4e2dd1a2ab"),
    ("cert_n29_L557.json", WB + "cert_n29_L557.json", "42b122994c40c02ceef410a314baf7d1c888ceb89e7a7b31af3c37de937654e0"),
    ("cert_n39_L650.json", WB + "cert_n39_L650.json", "dac136df114dd7a25cadb401bad63843dcf5be48d84af3cb9efd4d0ecccf66f3"),
]


def py(*args):
    subprocess.run([sys.executable, *args], check=True)


def main():
    out = sys.argv[1]; dl = os.path.join(out, "download"); orig = os.path.join(out, "orig")
    os.makedirs(dl, exist_ok=True); os.makedirs(orig, exist_ok=True)
    for name, url, sha in FILES:
        p = os.path.join(dl, name)
        if not os.path.exists(p):
            urllib.request.urlretrieve(url, p)
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if h != sha:
            sys.exit(f"SHA256 mismatch: {name} {h}")
        if name.endswith(".json"):  # wand125's JSON -> text (orig/wand_nXX_LYYY.txt)
            py(os.path.join(HERE, "json2txt.py"), p, os.path.join(orig, "wand_" + name[5:-5] + ".txt"))
        else:
            with open(p, "rb") as f, open(os.path.join(orig, name), "wb") as g:
                g.write(f.read())
    for f in sorted(os.listdir(orig)):
        py(os.path.join(HERE, "prep_cert.py"), os.path.join(orig, f), os.path.join(out, f))
    # our own negative controls
    py(os.path.join(HERE, "prep_cert.py"), os.path.join(orig, "s12_lower_3.931795_sparse.txt"),
       os.path.join(out, "mut_3.931795_sparse_scale.txt"), "--mutant", "scale", "999/1000")
    py(os.path.join(HERE, "prep_cert.py"), os.path.join(orig, "s12_lower_3.9686.txt"),
       os.path.join(out, "mut_3.9686_dropmax0.txt"), "--mutant", "dropmax", "0")
    # a mutant whose only violations are at θ > 0 (seed 1, 1,500 poses; about 3 minutes)
    py(os.path.join(HERE, "make_rot_mutant.py"), os.path.join(out, "s21_lower_4.9950.txt"),
       os.path.join(out, "mut_4.9950_rot.txt"), "--n", "1500")
    print("ready:", out)


if __name__ == "__main__":
    main()
