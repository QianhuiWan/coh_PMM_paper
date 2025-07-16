
import argparse
import re
from collections import defaultdict

def build_node_index(ins_node_map_file):
    node_to_ins = defaultdict(set)
    with open(ins_node_map_file) as f:
        for line in f:
            ins_id, nodes_str = line.strip().split('\t')
            for n in nodes_str.split(','):
                node_to_ins[n].add(ins_id)
    return node_to_ins

def main():
    parser = argparse.ArgumentParser(description="Extract detailed repeat-supporting read stats.")
    parser.add_argument("--gaf", required=True, help="Input GAF file")
    parser.add_argument("--node_map", required=True, help="INS-to-node mapping file")
    parser.add_argument("--out", required=True, help="Output TSV/txt file")
    parser.add_argument("--sample", required=True, help="Sample name for output")
    args = parser.parse_args()

    node_to_ins = build_node_index(args.node_map)

    with open(args.gaf) as f, open(args.out, "w") as out:
        out.write("sample\tins_id\tread_id\ttotal_node_hits\tunique_node_hits\tnodes_hit\n")
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            path = fields[5]
            nodes = [s[1:] for s in re.findall(r'[<>]\d+', path)]

            ins_hits = defaultdict(list)
            for node in nodes:
                for ins_id in node_to_ins.get(node, []):
                    ins_hits[ins_id].append(node)

            for ins_id, hit_nodes in ins_hits.items():
                total_hits = len(hit_nodes)
                unique_hits = len(set(hit_nodes))
                node_list = ','.join(hit_nodes)  # ordered, with duplicates
                out.write(f"{args.sample}\t{ins_id}\t{read_id}\t{total_hits}\t{unique_hits}\t{node_list}\n")

if __name__ == "__main__":
    main()


