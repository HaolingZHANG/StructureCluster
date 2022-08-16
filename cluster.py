from numpy import array, zeros, linalg, mean, sum

from monitor import Monitor


class Motif(object):

    def __init__(self, segment, structures):
        """
        Initialize the motif.

        :param segment: short amino acid sequence.
        :type segment: str

        :param structures: Referenced 3D structures.
        :type structures: list
        """
        self.segment = segment
        self.structures = structures


def collect_motifs(sequence, motif_length):
    """
    Collect motifs (by k-mer method) from all the amino acid sequence.

    :param sequence: amino acid sequence.
    :type sequence: list

    :param motif_length: length of clipping motif.
    :type motif_length: int

    :return: motif set generated in the protein database.
    :rtype: list

    .. note::
        reference: Phillip E. C. Compeau, Pavel A. Pevzner, and Glenn Tesler (2011) Nat. Biotechnol.
    """
    data, monitor = [], Monitor()
    if len(sequence) > motif_length + 1:
        for location in range(len(sequence) - motif_length + 1):
            data.append((sequence[location: location + motif_length], location))

    return data


def cluster_motifs(method, locations, position_groups, resolutions, verbose=False):
    """
    Cluster motif structures based on depth first search as connected components of undirected graphs.

    :param method: comparator for clustering structures.
    :type method: cluster.comparator.Similarity

    :param locations: locations in the protein database.
    :type locations: list

    :param position_groups: structures corresponded in each location in the protein database.
    :type position_groups: list

    :param resolutions: resolution of each corresponding structure in the protein database.
    :type resolutions: list

    :param verbose: need display the monitor.
    :type verbose: bool

    :return: cluster group of the token.
    :rtype: list

    .. note::
        reference: John Hopcroft, Robert Tarjan (1973) Commun. ACM
    """
    monitor = Monitor()

    if len(position_groups) == 1:
        return [[locations[0]]], position_groups[0], zeros(shape=(1, 1))

    if verbose:
        print("Cluster " + str(len(position_groups)) + " structures.")

    index_clusters, current = [], 0
    for candidate in range(len(locations)):
        related_indices = []
        for index, cluster in enumerate(index_clusters):
            for defender in cluster:
                result = method(structure_1=position_groups[candidate], resolution_1=resolutions[candidate],
                                structure_2=position_groups[defender], resolution_2=resolutions[defender])
                if result["similarity"]:
                    related_indices.append(index)
                    break
        if len(related_indices) == 0:
            index_clusters.append([candidate])
        elif len(related_indices) == 1:
            index_clusters[related_indices[0]].append(candidate)
        else:
            temp_clusters = []
            for index, cluster in enumerate(index_clusters):
                if (index not in related_indices) or (index == related_indices[0]):
                    temp_clusters.append(cluster)
                else:
                    temp_clusters[related_indices[0]] += cluster
            index_clusters = temp_clusters

        if verbose:
            monitor.output(current + 1, len(locations), extra={"C": len(index_clusters)})
            current += 1

    index_clusters = sorted(index_clusters, key=lambda c: len(c), reverse=True)

    if verbose:
        print("Generate reference structures.")

    clusters, reference_structures, reference_resolutions, inner_distances = [], [], [], []
    for index, index_cluster in enumerate(index_clusters):
        clusters.append([locations[location_index] for location_index in index_cluster])
        if len(index_cluster) > 1:
            distances = [0.0]
            structures = [position_groups[index_cluster[0]]]
            resolution = resolutions[index_cluster[0]]
            for location_index in index_cluster[1:]:
                result = method.calculate(structure_1=position_groups[location_index],
                                          resolution_1=resolutions[location_index],
                                          structure_2=position_groups[index_cluster[0]],
                                          resolution_2=resolutions[index_cluster[0]])
                structure = result["structures"][0]
                distance = linalg.norm(position_groups[index_cluster[0]] - structure, ord=2)
                distance /= len(position_groups[index_cluster[0]])
                distances.append(distance)
                structures.append(result["structures"][0])
                resolution += resolutions[location_index]
            structure = mean(array(structures), axis=0)
            resolution /= len(index_cluster)
            inner_distances.append(mean(array(distances)))
        else:
            structure = array(position_groups[index_cluster[0]])
            resolution = resolutions[index_cluster[0]]
            inner_distances.append(0.0)
        reference_structures.append(structure)
        reference_resolutions.append(resolution)

        if verbose:
            monitor.output(index + 1, len(index_clusters))

    if verbose:
        print("Generate atomic distances.")

    count = len(index_clusters)
    distances = zeros(shape=(count, count))
    current, total = 0, int((len(index_clusters) ** 2 - len(index_clusters)) / 2)
    for index_1 in range(count):
        distances[index_1][index_1] = inner_distances[index_1]
        for index_2 in range(index_1 + 1, count):
            distance = linalg.norm(reference_structures[index_1] - reference_structures[index_2], ord=2)
            distance /= len(reference_structures[index_1])
            distances[index_1][index_2] = distance
            distances[index_2][index_1] = distance

            if verbose:
                monitor.output(current + 1, total)
                current += 1

    return clusters, array(reference_structures), distances


def calculate_concentration(clusters):
    """
    Calculate the structure cluster concentration of a type of motif (1 minus Simpson's diversity index).
    |cite| Edward H. Simpson (1949) Nature

    :param clusters: structure clusters of a token.
    :type clusters: list

    :return: concentration index.
    :rtype: float
    """
    if len(clusters) == 1:
        return 1.0

    numbers = [len(cluster) for cluster in clusters]

    return sum([number * (number - 1) for number in numbers]) / (sum(numbers) * (sum(numbers) - 1))
