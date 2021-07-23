import argparse
import sys
import gzip


def read_metrics_file(file_name):
    metric_data = {}
    with gzip.open(file_name, "rt") if file_name.endswith(".gz") else open(
        file_name, "r"
    ) as fin:
        header = fin.readline()
        metric_data["header"] = header.split(",")[1:]
        for line in fin:
            fields = [x.strip() for x in line.split(",")]
            metric_data[fields[0]] = fields[1:]
    return metric_data


def compare_metric_files(ref_metrics, query_metrics):

    if ref_metrics["header"] != query_metrics["header"]:
        print("Header mismatch!")
        sys.exit(0)

    # print('Header', ref_metrics['header'])
    column_map_ref = {i: x.strip() for i, x in enumerate(ref_metrics["header"])}
    column_map_query = {i: x.strip() for i, x in enumerate(query_metrics["header"])}

    tol = 0.0000001
    nerrors = 0
    for key, fields in query_metrics.items():
        if key not in ref_metrics:
            print("Key {} missing in reference".format(key))
        if key != "header":
            # print()
            # print(key)
            # print('query: ', query_metrics[key])
            # print('ref  : ', ref_metrics[key])
            errors = False
            for i, (_x, _y) in enumerate(zip(query_metrics[key], ref_metrics[key])):
                if _x == "nan" or _y == "nan":
                    if _x != _y:
                        print(
                            "\tMismatch: key {}  query ({},  {}) and ref({}, {}) as  {}!={}".format(
                                key,
                                i,
                                column_map_query[i],
                                i,
                                column_map_query[i],
                                _x,
                                _y,
                            )
                        )
                        errors = True
                else:
                    x = float(_x)
                    y = float(_y)
                    if i not in [13, 14, 15, 16] and abs(x - y) > tol:
                        print(
                            "\tMismatch: key {}  query ({},  {}) and ref({}, {}) as  {}!={}".format(
                                key,
                                i,
                                column_map_query[i],
                                i,
                                column_map_query[i],
                                x,
                                y,
                            )
                        )
                        errors = True

            nerrors += int(errors)
            if nerrors > 0:
                sys.exit(0)


def main(args):
    # first create a list of sorted, and simplified sorted files
    ref_metrics = read_metrics_file(args.ref_file)
    query_metrics = read_metrics_file(args.query_file)
    i = 0

    print("Query against ref")
    compare_metric_files(ref_metrics, query_metrics)

    print("Ref against query")
    compare_metric_files(query_metrics, ref_metrics)


if __name__ == "__main__":
    description = """Compare the metric file """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--ref-file", dest="ref_file", required=True, help="Reference file of metrics"
    )
    parser.add_argument(
        "--query-file", dest="query_file", required=True, help="Query file of metrics"
    )
    args = parser.parse_args()

    main(args)
