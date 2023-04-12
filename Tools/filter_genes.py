
import csv


def find_genes(content: tuple, number: float, threshold: float, reverse: bool):
    """
    Finds up to the first `n` genes (there may not be that many) with 
    adj. p-value below `threshold`.

    The genes are ordered according to  its logFC and the `reverse` parameter.
    If true, they over-expressed genes are found. If False, from lowest 
    to highest under-expressed genes are found.
    """

    lines_with_na = 0
    lines_with_larger_adj_p_value = 0

    if reverse:
        filtered_content = tuple(
            filter(lambda gene: float(gene[4]) > 0, content[1:]))
    else:
        filtered_content = tuple(
            filter(lambda gene: float(gene[4]) < 0, content[1:]))

    sorted_content = sorted(
        filtered_content, key=lambda item: float(item[4]), reverse=reverse)

    expressed_genes = []

    for idx, line in enumerate(sorted_content):

        if idx == (number + lines_with_na + lines_with_larger_adj_p_value):
            break

        if 'NA' in line:
            lines_with_na += 1

            continue

        if float(line[8]) > threshold:
            lines_with_larger_adj_p_value += 1

            continue

        expressed_genes.append(line)

    return expressed_genes


def main(path: str, number: float, threshold: float) -> None:

    with open(path, 'r') as file:
        content = tuple(csv.reader(file, delimiter=','))

        over_expressed_genes = find_genes(content,
                                          number, threshold, reverse=True)
        under_expressed_genes = find_genes(content,
                                           number, threshold, reverse=False)

        with open('over_expressed.csv', 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(content[0])
            writer.writerows(over_expressed_genes)

        with open('under_expressed.csv', 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(content[0])
            writer.writerows(under_expressed_genes)


if __name__ == '__main__':

    import argparse

    ap = argparse.ArgumentParser()

    ap.add_argument('-p', '--path', required=True,
                    help='Path of microgravity.csv')
    ap.add_argument('-n', '--number', required=True,
                    help='Number of genes to save.')
    ap.add_argument('-t', '--threshold', required=True,
                    help='Highest adjusted p-value')

    args = vars(ap.parse_args())

    try:

        args['number'] = float(args['number'])
        args['threshold'] = float(args['threshold'])

        number = args['number']
        threshold = args['threshold']

    except ValueError:
        import sys
        print('Args are not valid')
        sys.exit(1)

    main(**args)
