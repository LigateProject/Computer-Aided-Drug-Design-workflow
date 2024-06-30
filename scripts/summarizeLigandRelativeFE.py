#!/usr/bin/python3
"""
This is a script to rank ligands based on pairwise relative FE differences.

The script assumes an input file, present in the directory where the script
is run, called ligandPairs.json. The file contains ligand pairs forming a
transformation sequence.

The script also assumes that there are directories (where the script is
executed) called "edge_<ligand_1_in_pair>_<ligand_2_in_pair>" each containing
a file RBFE.dat.

The script prints a ranked list of ligands and their relative free energies.
It also writes a JSON file with the results.

The simulations use a chain of pairwise transformations and yield their relative binding
free energies (RBFEs). There is also a set of "backup" pairs that can be used
if there is a break in the pairwise chain, e.g., if a simulation could not run
or finish. For each break in the chain the backup pair with the lowest uncertainty
is used to bridge the gap. If not all gaps can be filled, there will be separate
groups in the output, each listing their own internal relative free energies.

When all gaps are filled (that can be filled), there may still be backup pairs
that have not been used to fill any gaps. In that case, the output will present
the relative free energy difference of that backup pair using the chain of ligands
compared to if using the backup pair. That can give an estimation of how reliable
the results are.

Author: Magnus Lundborg
License: Apache License, version 2.0
"""

import copy
import json
import os

## Use the official uncertainties package if it is installed,
## otherwise fall back to the basic internalUFloat package.
#try:
    #from uncertainties import ufloat
    #USING_INTERNAL_UFLOAT = False
#except:
    #from internalUFloat import internalUFloat as ufloat
    #USING_INTERNAL_UFLOAT = True

# For development and testing force using the internalUFloat package.
# Does not keep track of correlated errors and may therefore overestimate
# the uncertainty.
from internal_ufloat import internal_ufloat as ufloat
USING_INTERNAL_UFLOAT = True



def extractLigandNameFromPath(path):
    """
    Extract the ligand name from a path, assuming a format similar to
    ./<ligand name>/*
    or
    <ligand name>/*
    """

    strippedPath = path.strip('./')
    parts = strippedPath.split('/')

    return parts[0]

def flattenListOfListByConcatenation(listOfLists):
    """
    Returns a list containing the flattened listOfLists, e.g.
    Input: listOfLists = [ [2, 4, 6, 8],
                           [1, 3, 5, 7],
                           [1, 2, 3, 4] ]
    Output: [2, 4, 6, 8, 1, 3, 5, 7, 1, 2, 3, 4]
    """
    flattenedList = []
    for entry in listOfLists:
        flattenedList += entry
    return flattenedList

def getRelDiffOfPairFromFile(ligand1Name, ligand2Name,
                             directoryPrefix = 'edge', fileName = 'RBFE.dat'):
    """
    Get the relative binding free energy difference in a ligand pair.
    The RBFEs are read from the file in each ligand pair directory (with a directory prefix).
    The values are returned as internalUFloat if there is an uncertainty value, otherwise float.
    """

    dirName = '_'.join([directoryPrefix, ligand1Name, ligand2Name])
    filePath = os.path.join(dirName, fileName)
    #print(filePath)
    result = None
    try:
        with open(filePath) as f:
            for line in f:
                parts = line.split()
                # Accept either of
                # <relative difference> (returns a float)
                # <relative difference> <uncertainty> (returns an internalUFloat)
                # <relative difference> +- <uncertainty> (returns an internalUFloat)
                if len(parts) == 1:
                    try:
                        result = float(parts[0])
                        break
                    except ValueError:
                        pass
                elif len(parts) == 2 or len(parts) == 3 and parts[1] == '+-':
                    try:
                        # The arguments must be converted to float before making the internalUFloat.
                        diff = float(parts[0])
                        uncertainty = float(parts[-1])
                        result = ufloat(diff, uncertainty)
                        break
                    except ValueError:
                        pass
            else:
                print(f'Warning: Cannot find relative difference between pair '
                    f'{ligand1Name} and {ligand2Name}.')

    except FileNotFoundError:
        print(f'Cannot find {filePath}. Continuing without RBFE data for '
              f'{ligand1Name} and {ligand2Name}.')

    return result

def makeLigandPairListsWithRelativeDifference(inputLigandPairList):
    """
    Return a list of continuous ligand pair lists, including their pairwise RBFEs,
    and one with backup ligand pairs.

    If there are non-existing RBFEs there will be gaps in the ligand pair list. These
    are divided into separate lists in this function. If there are pairwise RBFEs
    for all ligands, the returned list will only contain one list of ligand pairs.
    """

    listOfLigandPairLists = []
    currentListOfLigandPairs = []
    backupLigandPairs = []
    # Iterate through the pairs of ligands. Keep the list of ligand pairs separated from a list of
    # potential backup ligand pairs. The backup ligand pairs are always placed last in the input
    # list.
    # The backup ligand pairs already exist in the list of ligand pairs.
    for inputPair in inputLigandPairList:
        newPair = {}
        newPair['Ligand_1'] = extractLigandNameFromPath(inputPair['Ligand_1'])
        newPair['Ligand_2'] = extractLigandNameFromPath(inputPair['Ligand_2'])
        relativeDifference = getRelDiffOfPairFromFile(newPair['Ligand_1'], newPair['Ligand_2'])
        newPair['pairwise_difference'] = relativeDifference
        newPair['relative_fe_difference'] = ufloat(0,0) # This will be used later

        # Here we need to use a flattened list of all lists to check all added pairs.
        flattenedListOfLigandPairs = (flattenListOfListByConcatenation(listOfLigandPairLists) +
                                      currentListOfLigandPairs)

        # Check if either of the ligands in the pair is already present as first ligand in the
        # list of pairs. In that case this is a backup pair.
        if ((len(flattenedListOfLigandPairs) > 0) and
                (next((item for item in flattenedListOfLigandPairs if item['Ligand_1'] ==
                       newPair['Ligand_1']), None) is not None or
                 next((item for item in flattenedListOfLigandPairs if item['Ligand_1'] ==
                       newPair['Ligand_2']), None) is not None)):
            # Do not append backup pairs that do not have any RBFE results.
            if relativeDifference is not None:
                backupLigandPairs.append(newPair)
        else:
            # Start a new list if there is a break in the chain (i.e., no RBFE value in a pair)
            if relativeDifference is None:
                # If there were no entries in the list beforehand, add this entry before creating a
                # new list to avoid losing any ligand information from this missing pair.
                # It is not important that all pairs with no pairwise RBFEs are stored, but all
                # ligands must be recorded.
                if not currentListOfLigandPairs:
                    currentListOfLigandPairs.append(newPair)
                listOfLigandPairLists.append(currentListOfLigandPairs)
                currentListOfLigandPairs = []
                continue
            currentListOfLigandPairs.append(newPair)

    if currentListOfLigandPairs:
        listOfLigandPairLists.append(currentListOfLigandPairs)


    return (listOfLigandPairLists, backupLigandPairs)

def findBestBackupPair(firstList, secondList, backupLigandPairs):
    """
    Find the best backup pair to bridge a gap between the two lists.
    """

    bestPair = None
    for pair in backupLigandPairs:
        # Check if both ligands of the ligand pair are present as one of two ligands in a pair in both lists.
        if (((next((item for item in firstList if item['Ligand_1'] == pair['Ligand_1']), None) is not None) or
             (next((item for item in firstList if item['Ligand_2'] == pair['Ligand_1']), None) is not None)) and
            ((next((item for item in secondList if item['Ligand_1'] == pair['Ligand_2']), None) is not None) or
             (next((item for item in secondList if item['Ligand_2'] == pair['Ligand_2']), None) is not None))):
            if bestPair is None or pair['pairwise_difference'].std_dev < bestPair['pairwise_difference'].std_dev:
                bestPair = pair
        elif (((next((item for item in secondList if item['Ligand_1'] == pair['Ligand_1']), None) is not None) or
               (next((item for item in secondList if item['Ligand_2'] == pair['Ligand_1']), None) is not None)) and
              ((next((item for item in firstList if item['Ligand_1'] == pair['Ligand_2']), None) is not None) or
               (next((item for item in firstList if item['Ligand_2'] == pair['Ligand_2']), None) is not None))):
            if bestPair is None or pair['pairwise_difference'].std_dev < bestPair['pairwise_difference'].std_dev:
                bestPair = pair

    return bestPair

def findPairlistBackupPairConnections(listOfLigandPairLists, backupLigandPairs):
    """
    Try to find backup pairs to bridge gaps between the separated lists in
    listOfLigandPairLists.

    A gap will be filled by using the backup pair with the lowest uncertainty,
    if there is more than one matching backup pair.

    Returns a list of backup pairs matching the gaps in the listOfLigandPairLists.
    If there was no backup pair found for a gap it will be None
    and
    a list of connection tuples where (0, 1) means that the backup pair of the same
    index was used to connect the (previously disconnected lists) 0 and 1 in
    listOfLigandPairLists.
    """

    usedBackupLigandPairs = []
    connections = []

    for i, outerList in enumerate(listOfLigandPairLists[:-1]):
        for j, matchingList in enumerate(listOfLigandPairLists[(i + 1):]):
            bestPair = findBestBackupPair(outerList, matchingList, backupLigandPairs)

            if bestPair is not None and bestPair not in usedBackupLigandPairs:
                # Check direction of the connection.
                if next((item for item in outerList if item['Ligand_1'] == bestPair['Ligand_1']), None) is not None:
                    connections.append((i, j + 1))
                else:
                    connections.append((j + 1, i))
                # Keek track of which lists have been connected using a backup pair.
                usedBackupLigandPairs.append(bestPair)
            else:
                usedBackupLigandPairs.append(None)
                connections.append(None)

    return usedBackupLigandPairs, connections


def relativeDifferenceBetweenPairs(pairList, ligand1Name, ligand2Name):
    """
    Calculate and return the relative FE difference between two ligands (by name)
    """

    reverseOrder = False

    ligand1 = next((item for item in pairList if item['Ligand_1'] == ligand1Name), None)
    if ligand1 is not None:
        ligand2 = next((innerItem for innerItem in pairList if innerItem['Ligand_2'] ==
                        ligand2Name), None)
    else:
        ligand1 = next((item for item in pairList if item['Ligand_1'] == ligand2Name), None)
        if ligand1 is not None:
            ligand2 = next((innerItem for innerItem in pairList if innerItem['Ligand_2'] ==
                            ligand1Name), None)
            reverseOrder = True

    if ligand1 is None or ligand2 is None:
        print(f'Cannot find endstate connections between {ligand1Name} and {ligand2Name}')
        return None

    ligandPair = ligand1
    if USING_INTERNAL_UFLOAT:
        relativeDifference = ufloat(ligandPair['pairwise_difference'])
    else:
        relativeDifference = ligandPair['pairwise_difference']
    while ligandPair != ligand2 and ligandPair is not None:
        ligandPair = next((item for item in pairList if item['Ligand_1'] ==
                           ligandPair['Ligand_2']), None)
        if ligandPair is not None:
            try:
                relativeDifference += ligandPair['pairwise_difference']
            except (TypeError, ValueError):
                return None

    if reverseOrder:
        relativeDifference *= -1

    return relativeDifference

def applyRelativeDifferenceToPairLists(pairLists, backupPairs, connections):
    """
    Iterate through all connections. Find the connected pair lists. Apply a correction so that the relative difference
    between the lists (connected by a backup pair) matches the backup pair.
    """

    for i, backupPair in enumerate(backupPairs):
        connection = connections[i]
        if backupPair is None or connection is None:
            if (backupPair is None and connection is not None) or (connection is None and backupPair is not None):
                print('Warning! Mismatch in backup pairs and connections')
            continue
        fromList = connection[0]
        toList = connection[1]
        fromLigandName = backupPair['Ligand_1']
        toLigandName = backupPair['Ligand_2']
        relativeDifferenceInConnection = backupPair['pairwise_difference']

        # The pairwise relative differences in the pair lists are stored in ligand_2.
        fromLigand = next((item for item in pairLists[fromList] if item['Ligand_2'] == fromLigandName), None)
        toLigand = next((item for item in pairLists[toList] if item['Ligand_2'] == toLigandName), None)

        if fromLigand and toLigand:
            appliedRelativeDifference = (fromLigand['relative_fe_difference'] +
                                        relativeDifferenceInConnection -
                                        toLigand['relative_fe_difference'])

            if toList:
                for ligand in pairLists[toList]:
                    ligand['relative_fe_difference'] += appliedRelativeDifference


def calculateRelativeDifferenceToFirstLigand(ligandPairList):
    """
    For ligand in a list calculate the relative difference to the first ligand in the list.
    The relative difference is stored in the pair dictionary as 'relative_fe_difference'.
    """

    firstPair = ligandPairList[0]
    for item in ligandPairList:
        relativeDifferenceToFirstLigand = relativeDifferenceBetweenPairs(ligandPairList,
                                                                         firstPair['Ligand_1'],
                                                                         item['Ligand_2'])
        item['relative_fe_difference'] = relativeDifferenceToFirstLigand

    # The list above is created based on "accumulated" pairwise free energy differences compared
    # to the first entry in the input.
    # Add a dummy entry to make sure that the first ligand is shown in the sorted output list.
    ligandPairList.append({'Ligand_1' : firstPair['Ligand_1'],
                           'Ligand_2' : firstPair['Ligand_1'],
                           'pairwise_difference' : ufloat(0, 0),
                           'relative_fe_difference' : ufloat(0, 0)})

def mergeLigandPairListsBasedOnConnections(listOfLigandPairLists, backupPairs, connections):
    """
    Merge the pair lists that are connected by backup pairs

    Returns true if any lists were merged.
    """

    # Iterate through the connections in reverse order with an enumerated list

    hasMergedAnyLists = False
    indicesToRemove = []

    enumeratedBackupPairs = tuple(enumerate(backupPairs))
    for i, backupPair in reversed(enumeratedBackupPairs):
        connection = connections[i]
        if backupPair is None or connection is None:
            if (backupPair is None and connection is not None) or (connection is None and backupPair is not None):
                print('Warning! Mismatch in backup pairs and connections.')
            continue

        sourceListIndex = max(connection)
        destinationListIndex = min(connection)

        if sourceListIndex in indicesToRemove or destinationListIndex in indicesToRemove:
            print('Error! Trying to merge lists that are no longer valid. Aborting list merging. Results may not be valid.')
            return hasMergedAnyLists

        sourceList = listOfLigandPairLists[sourceListIndex]
        destinationList = listOfLigandPairLists[destinationListIndex]

        #print(f'Merging {sourceList} into {destinationList}')

        hasMergedAnyLists = True
        indicesToRemove.append(sourceListIndex)
        destinationList += sourceList

    # Prune the list of lists after all merges to keep the indices valid.
    for i in indicesToRemove:
        del listOfLigandPairLists[i]

    return hasMergedAnyLists

def printComparisonToUnusedBackupPairs(listOfLigandLists, backupPairs):
    """
    For all backupPairs within a fully connected list of ligands,
    print a comparison of the relative BFE of the result in the list vs
    the BFE of the backupPair.
    """

    if not backupPairs:
        return

    print('Comparing results to backup pairs:')
    for backupPair in backupPairs:
        fromLigandName = backupPair['Ligand_1']
        toLigandName = backupPair['Ligand_2']
        for listOfLigands in listOfLigandLists:
            fromLigand = next((item for item in listOfLigands if item['ligand'] == fromLigandName), None)
            toLigand = next((item for item in listOfLigands if item['ligand'] == toLigandName), None)
            if fromLigand is not None and toLigand is not None:
                relativeBfe = toLigand['relative_fe_difference'] - fromLigand['relative_fe_difference']
                print(f'{fromLigandName} to {toLigandName}: Relative BFE from list: {relativeBfe} vs from backup pair: {backupPair["pairwise_difference"]}')

    print('\n')


def main():
    """
    The main function of the script.

    - Read a JSON file ('ligandPairs.json') with input pairs.
    - Create a list of pairs including their pairwise free energy differences,
      if there are gaps make multiple lists of unbroken chains.
    - Calculate the "accumulated" free energy difference of each ligand to the
      first ligand in the input list.
    - Fill gaps using backup pairs, if possible. Gaps are bridged using the backup
      pair, that connects two lists, with the lowest uncertainty.
    - Sort the list(s) and recalibrate RBFE 0 for the highest ranked ligand, in each
      disconnected list.
    - Print the relative free energies of all disconnected ligand lists.
    - Make a comparison to the free energies from backup pairs that have not been used
      to bridge gaps.
    - Output a JSON file ('ligandRanking.json')
    """

    with open('ligandPairs.json') as inputFile:
        inputLigandPairs = json.load(inputFile)

    (listOfLigandPairLists, backupLigandPairs) = makeLigandPairListsWithRelativeDifference(inputLigandPairs)

    print(f'Backup pairs: {backupLigandPairs}')

    (usedBackupPairs, gapConnections) = findPairlistBackupPairConnections(listOfLigandPairLists, backupLigandPairs)
    print(f'Backup pairs used to make connections between disconnected pair list segments: {usedBackupPairs}\n')

    for i, ligandPairList in enumerate(listOfLigandPairLists):
        if any(ligandPair['pairwise_difference'] is None for ligandPair in ligandPairList):
            print(f'Warning: Pair list contains a pair with an invalid free energy difference. Ignoring this block: {listOfLigandPairLists[i]}.')
            del listOfLigandPairLists[i]

    for ligandPairList in listOfLigandPairLists:
        calculateRelativeDifferenceToFirstLigand(ligandPairList)

    applyRelativeDifferenceToPairLists(listOfLigandPairLists, usedBackupPairs, gapConnections)

    listOfSortedLigandLists = []
    for ligandPairList in listOfLigandPairLists:
        sortedLigandPairs = sorted(ligandPairList, key = lambda x: x['relative_fe_difference'])

        # From the ligand pairs, make a clean list of ligands and their relative_fe_difference
        sortedLigands = []
        for pair in sortedLigandPairs:
            sortedLigands.append({'ligand' : pair['Ligand_2'],
                                'relative_fe_difference' :
                                    pair["relative_fe_difference"]})

        listOfSortedLigandLists.append(sortedLigands)

    needsResorting = mergeLigandPairListsBasedOnConnections(listOfSortedLigandLists, usedBackupPairs, gapConnections)
    if needsResorting:
        for i, ligandList in enumerate(listOfSortedLigandLists):
            listOfSortedLigandLists[i] = sorted(ligandList, key = lambda x: x['relative_fe_difference'])

    # Adjust the free energies so that 0 is the lowest in each disconnected list
    for sortedList in listOfSortedLigandLists:
        lowestRelDiff = sortedList[0]['relative_fe_difference']
        nominalDiffValue = float(lowestRelDiff.nominal_value)
        for item in sortedList:
            item['relative_fe_difference'] -= nominalDiffValue

    # Print the output
    for i, sortedList in enumerate(listOfSortedLigandLists):
        print(f'Ligand chain {(i + 1)}:')
        for item in sortedList:
            print(f'Ligand: {item["ligand"]}\t:\t{item["relative_fe_difference"]}')
        print()

    # Print a comparison to backup pairs that have not been used to merge lists. This can help evaluate if the results would
    # be the same if another path is taken.
    printComparisonToUnusedBackupPairs(listOfSortedLigandLists, [item for item in backupLigandPairs if item not in usedBackupPairs])

    # Convert ufloat values to strings for JSON output
    for sortedList in listOfSortedLigandLists:
        for item in sortedList:
            item['relative_fe_difference'] = f'{item["relative_fe_difference"]}'

    # Write to a JSON output file
    outputFileName = 'ligandRanking.json'
    with open(outputFileName, 'w') as outputFile:
        if len(listOfSortedLigandLists) > 1:
            json.dump(listOfSortedLigandLists, outputFile, indent=2)
        else:
            json.dump(listOfSortedLigandLists[0], outputFile, indent=2)

    print(f'Results written to {outputFileName}\n')

if __name__ == '__main__':
    main()
