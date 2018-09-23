"""CSCA08 Assignment 2, Fall 2017
 I hereby agree that the work contained herein is solely my work and that I
 have not received any external help from my peers, nor have I used any
 resources not directly supplied by the course in order to complete this
 assignment. I have not looked at anyone else's solution, and no one has
 looked at mine. I understand that by adding my name to this file, I am
 making a formal declaration, and any subsequent discovery of plagiarism
 or other academic misconduct could result in a charge of perjury in
 addition to other charges under the academic code of conduct of the
 University of Toronto Scarborough Campus
 Name: Chou Lu Wei
 UtorID: Choulu1
 Student Number: 1004295693
 Date: 2017/12/03
"""
# define the how many chromosome pair a human has. we start counting from 0
# to 22 inclusive
human_chromosome = 22
# define the location of the human chromosome sex pair which is 22nd chromsome
h_chromo_sex = 22
# list keeping track of all human id no two id can be the same
human_id_list = []
# define the the memory nucleotides
m_nucleo_sym = '0123456789'
# define the left maternal and right maternal symbol
l_mat_sym = 'LM'
r_mat_sym = 'RM'
# define symbol used to represent male and female
male_symbol = 'M'
female_symbol = 'F'
# define variable to keep track of number of child created used for creating
# unique id for each child
child_born = 0


class PairNotFoundError(Exception):
    ''' Raised when pair was not found for a given position number
    '''


class ChromosomeNotFoundError(Exception):
    ''' Raised when Chromosome Pair was not found for a given pair number
    '''


class MarkerNotFoundError(Exception):
    ''' Raised when Marker is not found in a Chromosome Set
    '''


class NullPairError(Exception):
    ''' Raised when user initialize marker but does not set towards anything
    and the user tries to access the empty marker
    '''


class SameIdError(Exception):
    ''' Raised when two things have the same id
    '''


class SexlessBabyError(Exception):
    ''' Raised when Binder given does not have a sex when famale procreate
    '''


class InsuffientGeneError(Exception):
    ''' Raised when females procreate and either Male or Female in
    procreate method does not have gene specified at a index in the binder
    '''


class ChromosomePair():
    ''' A object representing a chromosome pair. The majority human DNA are
    identical and the chromosome pair class will only records difference in
    genes between individuals
    '''
    def __init__(self):
        ''' (ChromosomePair) -> NoneType
        Initializes a ChromosomePair object
        '''
        # initialize dictionary which stores maps the location to the pair
        self._pos_to_pair = dict()

    def set_by_pos(self, pos, pair):
        ''' (ChromosomePair, int, str) -> NoneType
        REQ: pos >= 0
        REQ: len(pair) == 2
        Returns none given a position integer and string pair, stores the
        position and the pair in the chromosome class
        '''
        # store the position in our object
        self._pos_to_pair[pos] = pair

    def get_by_pos(self, pos):
        ''' (ChromosomePair, int) -> str
        REQ: pos >= 0
        RAISE: PairNotFoundError
        Returns a pair given an integer position.
        - Raises PairNotFoundError if no pair was found given the position
        '''
        # try to return a pair given a position by accessing our dictionary
        try:
            return self._pos_to_pair[pos]
        # if we can't access position raise pos error
        except KeyError:
            raise PairNotFoundError("Pair was not found for at position " +
                                    str(pos) + " within the ChromosomePair")

    def get_all_pos(self):
        '''(ChromosomePair) -> set
        Returns all the positions within Chromosome Pair as a set like object
        '''
        # use key method for dictionary to return all positions
        return self._pos_to_pair.keys()


class ChromosomeSet():
    ''' Class representing multiple chromosome pairs
    '''
    def __init__(self, total_pair, sex_chromo_pos):
        ''' (ChromosomeSet, int, int) -> NoneType
        REQ: total_pair and sex_chromo_pos >= 0
        Intialize a ChromosomeSet object given an ChromosomeSet and an
        integer which defines how many pairs of Chromosome is in set and
        the position of the sex chromosome
        - Raise InvalidInputError if total pair or sex chromosome is negative
        number
        '''
        # initialize the num to chromosome dictionary
        self._num_to_chromosome = dict()
        # initialize the marker to pair location dictionary
        self._marker_to_pair = dict()
        # initialize variable defining total pair
        self._total_pair = total_pair
        # initialize variable defining where the sex chromosome is
        self._sex_chromo = sex_chromo_pos
        # start a counter to loop from 0 to total pair + 1 we add 1 since we
        # want to include all values between and including 0 and total_pair
        for i in range(0, total_pair + 1):
            # set the i-th key of the dictionary to have the i-th chromosome
            i_chromosome = ChromosomePair()
            self._num_to_chromosome[i] = i_chromosome

    def get_total_pair(self):
        ''' (ChromosomeSet) -> int
        Return the total number of Chromosome Pair in set given ChromosomeSet
        Note: we do allow user to set total pair because setting total pairs
        would fundmentally change the structure of the data
        '''
        return self._total_pair

    def get_sex_chromo(self):
        ''' (ChromosomeSet) -> int
        Return sex chromosome position in set given ChromosomeSet
        '''
        return self._sex_chromo

    def set_sex_chromo(self, new_pos):
        ''' (ChromosomeSet, int) -> int
        REQ: 0 <= new_pos <= self.get_total_pair()
        sets the sex chromosome given a chromsome position
        '''
        self._sex_chromo = new_pos

    def get_chromosome(self, pair_num):
        ''' (ChromosomeSet, int) -> ChromosomePair
        RAISE: ChromosomeNotFoundError
        REQ: pair_num >= 0
        Return a ChromosomePair given a ChromosomeSet and an integer pair
        number.
        - Raises ChromosomeNotFoundError if no ChromosomePair was found
        in ChromosomeSet object
        - Raise InvalidInputError if pair_number is less than 0
        '''
        # try to return a ChromosomePair if no such ChromosomePair exist
        # raise ChromosomeNotFoundError
        try:
            return self._num_to_chromosome[pair_num]
        except KeyError:
            raise ChromosomeNotFoundError("Chromosome Pair was not found at" +
                                          "pair: " + str(pair_num) +
                                          " given location of index")

    def set_chromosome(self, pair_num, new_chromo):
        '''(ChromosomeSet, int, ChromosomePair) -> NoneType
        REQ: pair_num >= 0
        Sets the Chromosome Pair at pair number to be a New ChromosomePair
        - Raise InvalidInputError if pair_number is less than 0
        '''
        self._num_to_chromosome[pair_num] = new_chromo

    def set_by_pos(self, pair_num, pos, new_pair):
        '''(ChromosomeSet, int, int, str) -> NoneType
        RAISE: ChromosomeNotFoundError
        REQ: len(new_pair) == 2
        REQ: pair_num and pos >= 0
        Sets the at the pair number and position number to be a
        new pair.
        - Raises ChromosomeNotFoundError if Pair is not found at index
        '''
        # get the chromosome using get chromosome method
        chromosome_i = self.get_chromosome(pair_num)
        # once we have access to chromosome use set by pos in ChromosomePair
        # to access our pairs and set a new pair
        chromosome_i.set_by_pos(pos, new_pair)

    def get_by_pos(self, pair_num, pos):
        '''(ChromosomeSet, int, int) -> str
        REQ: pair_num and pos >= 0
        RAISE: ChromosomeNotFoundError
        RAISE: PairNotFoundError
        Returns a pair represented as string given a integer
        pair number and integer position number.
        - Raises ChromosomeNotFoundError if chromosome pair is not found.
        - Raises PairNotFoundError if pair is not found at index
        '''
        # try to access the chromosome at the pair number using get chromosome
        # method. Method will Raise ChromosomeNotFound if Chromosome is not
        # was not found at pair number
        chromosome_i = self.get_chromosome(pair_num)
        # access exact pair using get by pos method, method will
        # Raise PairNotFoundError if pair is not found
        pair_i = chromosome_i.get_by_pos(pos)
        return pair_i

    def set_marker(self, marker, pair_num, pos):
        '''(ChromosomeSet, str, int, int) -> NoneType
        REQ: pos >= 0
        RAISE: ChromosomeNotFoundError
        Sets a marker to point to the position of a pair given a pair number
        and position number.
        - Raises ChromosomeNotFoundError if pair number was not found in
        ChromosomeSet object
        '''
        # Try to access the chromosome pair at which the specified index
        # through the use of the get by pos method which can raise
        # PairNotFoundError or ChromosomeNotFoundError.
        try:
            self.get_by_pos(pair_num, pos)
        # if chromosome was not found then marker is pointing to an invalid
        # location and will method will raise ChromosomeNotFoundError
        # if pair was not found at location
        except PairNotFoundError:
            # except the error and do nothing user will be require to set it
            # later with set by marker function
            pass
        # if both chromosome and pair was found or only chromosome was found
        # set the marker to point the index
        self._marker_to_pair[marker] = (pair_num, pos)

    def _get_by_marker_index(self, marker):
        '''(ChromosomeSet, str) -> (int, int)
        RAISE:  MarkerNotFoundError
        Return a tuple of integer representing indexes given a marker which
        points towards a pair in ChromosomeSet.
        - Raises MarkerNotFoundError if no marker was found in the
        ChromosomeSet object
        '''
        # try to access the marker to find what it points to
        try:
            return self._marker_to_pair[marker]
        # if marker is not found in our dictionary raise MarkerNotFoundError
        except KeyError:
            raise MarkerNotFoundError('Marker: ' + str(marker) +
                                      ' was not found in a ChromosomeSet')

    def set_by_marker(self, marker, new_pair):
        ''' (ChromosomeSet, str, str) -> NoneType
        REQ: len(new_pair) == 2
        RAISE: MarkerNotFoundError
        Sets pair pointed at by Marker to be a new pair.
        - Raises MarkerNotFoundError if no marker was found in ChromosomeSet
        object
        '''
        # use get marker method to acquire the index of the position we want
        # method will Raise MarkerNotFoundError if marker does not exist
        # within our object
        (pair_num, pos) = self._get_by_marker_index(marker)
        # use set by pos method to change to new pair
        self.set_by_pos(pair_num, pos, new_pair)

    def get_by_marker(self, marker):
        ''' (ChromosomeSet, str) -> str
        RAISE: MarkerNotFoundError
        Returns a pair given the marker which the points towards the pair.
        - Raises MarkerNotFoundError if no marker was found in ChromosomeSet
        object.
        '''
        # use get marker method to acquire the index. get marker method will
        # throw MarkerNotFoundError if marker was not found
        (pair_num, pos) = self._get_by_marker_index(marker)
        # try to use get by pos method to acquire the return value.
        try:
            return_value = self.get_by_pos(pair_num, pos)
        # if we cannot find the pair
        except PairNotFoundError:
            # raise NullPairError
            raise NullPairError('Unable to retrieve pair because marker ' +
                                marker + ' was not specified to point' +
                                ' towards a value')
        return return_value


class Human(ChromosomeSet):
    ''' A class representing a human. In this case we represent a Human
    as a set of 23 chromosome pairs with an id.
    '''
    def __init__(self, human_id):
        ''' (Human, str) -> NoneType
        RAISE: SameIdError
        Intialize a human object
        - Raises SameIdError if ID of human object is the same as another
        previously defined object
        '''
        # define our global variable which is human_chromosome, h_chro
        # initialize parent class and set it to have a set of chromosome
        # with human_chromosome number of chromosome
        ChromosomeSet.__init__(self, human_chromosome, h_chromo_sex)
        # initialize human id and store it in the object
        self._human_id = human_id
        # check if id is in global human id list by iterating through global
        # human list
        if (human_id in human_id_list):
            # if it is raise Same Id Error
            raise SameIdError('ID: ' + human_id + ' already exist. No two' +
                              ' Humans can have the same ID')
        else:
            # otherwise just append it to the global list
            human_id_list.append(human_id)

    def get_id(self):
        ''' (Human) -> str
        Return Human Id given Human
        '''
        return self._human_id

    def set_id(self, new_id):
        '''(Human) -> None
        RAISE: SameIdError
        Change id of human
        - Raises SameIdError if new ID of human object is the same as another
        previously defined object
        '''
        # remove id from global list because we want a new id
        human_id_list.remove(self._human_id)
        # check if new id is in our global human list
        if (new_id in human_id_list):
            # if it is since we can't set two id to be the same we raise
            # SameIdError
            raise SameIdError('ID: ' + new_id + ' already exist. No two' +
                              ' Humans can have the same ID')
        else:
            # otherwise set a new id
            self._human_id = new_id
            # add that new id to our global list
            human_id_list.append(new_id)

    def set_by_pos(self, pair_num, pos, new_nucleo_pair):
        ''' (Human, int, int, str) -> NoneType
        REQ: new_nucleo_pair must only contain nucleotides by default
        these nucleotides are 'ATCG'
        REQ: pair_num and pos >= 0
        RAISE: ChromosomeNotFoundError
        Sets the nucleotides at the pair number and position number to be a
        new pair of nucleotides.
        - Raises ChromosomeNotFoundError if Chromosome Pair is not found
        '''
        # We overwrite parent method for new requirement
        ChromosomeSet.set_by_pos(self, pair_num, pos, new_nucleo_pair)

    def set_by_marker(self, marker, new_nucleo_pair):
        ''' (Human, str, str) -> NoneType
        REQ: new_nucleo_pair must only contain nucleotides by default
        these nucleotides are 'ATCG'
        RAISE: MarkerNotFoundError
        Sets nucleotide pointed at by Marker to be a new pair of nucleotides.
        - Raises MarkerNotFoundError if no marker was found in ChromosomeSet
        object
        '''
        ChromosomeSet.set_by_marker(self, marker, new_nucleo_pair)

    def test(self, query):
        ''' (Human, Query) -> bool
        Return a boolean indicating if a Query matched a Human's chromosome
        given a Query and a Human.
        '''
        # define our globals variables to be m_nucleo_sym which is what the
        # the program defines as memory nucleotides
        global m_nucleo_sym
        # define a dictionary for the memory nucleotides
        mem_to_nucleo = dict()
        # define if Human can compare with human chromosome which by default we
        # assume to be true
        can_compare = True
        # define the pair number we are on which is zero by default since
        # we start counting from 0
        pair_num = 0
        # define what our total chromosome is and where our sex chromosome is
        human_chromosome = self.get_total_pair()
        h_chromo_sex = self.get_sex_chromo()
        # start while loop which loops only when we can compare the Query
        # Values and when the pair we are comparing is less than the total
        # human chromosome
        while (can_compare and pair_num <= human_chromosome):
            # get each Chromosome Pair positions
            pos_set = query.get_chromosome(pair_num).get_all_pos()
            # loop through all the position in the Chromosome Pairs
            for pos in pos_set:
                # define the actual compare value to be empty string this is
                # the nucleotide with all the memory nucleotide translated
                q_compare = ''
                # define what is at index in query
                q_nucleo = query.get_by_pos(pair_num, pos)
                # try to access human nucleotide a pair number and pos
                try:
                    h_nucleo = self.get_by_pos(pair_num, pos)
                    # loop through the individual nucleotide pairs
                    for index in range(len(q_nucleo)):
                        # if a nucleotide is a memory nucloetide
                        # and it not in our memory nucleotide dictionary
                        if ((q_nucleo[index] in m_nucleo_sym) and
                                (q_nucleo[index] not in mem_to_nucleo)):
                            # add it to our dictionary with it pointing to
                            # the nucleotide from human
                            mem_to_nucleo[q_nucleo[index]] = h_nucleo[index]
                            # and then add it to our actual compare value
                            q_compare += h_nucleo[index]
                        # otherwise if number in the dictionary
                        elif (q_nucleo[index] in mem_to_nucleo):
                            # look up number and add it to our actual
                            # compare value
                            q_compare += mem_to_nucleo[q_nucleo[index]]
                        # else it must be a nucleotide
                        else:
                            # in which case we add to to our actual compare
                            # value
                            q_compare += q_nucleo[index]
                    # if our memory query does not equal to human nucleotide
                    if (q_compare != h_nucleo):
                        # then we can set can_compare to be false
                        can_compare = False
                # if we cannot find Pair except PairNotFoundError
                except (PairNotFoundError):
                    # if the pair is chromosome sex pair and the self is Male
                    if (pair_num == h_chromo_sex and type(self) is Male):
                        # then we set can compare to be false
                        can_compare = False
            # increment pair num to move to next chromosome
            pair_num += 1
        # return can_compare for the result of the whole query
        return can_compare


class Query(ChromosomeSet):
    ''' A class which represents a query. A query is a special Chromosome Set
    acts like a mask to check if certain patterns emerge from a set of Human
    chromosomes. For future reference this Query classes only deals with Human
    chromosomes it might be useful to rename this as Human_Query if other life
    forms are implemented.
    '''
    def __init__(self):
        ''' (Query) -> Nonetype
        Initializes a Query object
        '''
        # initialize how many chromosomes a human has and define pass in
        # where the sex chromsome is
        ChromosomeSet.__init__(self, human_chromosome, h_chromo_sex)

    def set_by_pos(self, pair_num, pos, query_pair):
        ''' (Query, int, int, str) -> NoneType
        REQ: pair_num and pos >= 0
        REQ: query_pair can only contain nucleotides or memory nucleotides
        by default those are '0123456789' or 'ATCG'
        RAISE: ChromosomeNotFoundError
        Sets nucleotide pointed at by Marker to be a new pair of nucleotides.
        nucleotide included can also be special memory nucleotide which are
        denoted as string numbers 0 to 9.
        - Raises ChromosomeNotFoundError if Chromosome Pair is not found
        '''
        # overwrite parent method for new requirements
        ChromosomeSet.set_by_pos(self, pair_num, pos, query_pair)

    def set_by_marker(self, marker, query_pair):
        ''' (Query, str, str) -> NoneType
        REQ: query_pair can only contain nucleotides or memory nucleotides
        by default those are '123456789' or 'ATCG'
        RAISE: MarkerNotFoundError
        Sets nucleotide pointed at by Marker to be a new pair of nucleotides.
        nucleotide included can also be special memory nucleotide which are
        denoted as string numbers 0 to 9.
        - Raises MarkerNotFoundError if no marker was found in ChromosomeSet
        object
        '''
        ChromosomeSet.set_by_marker(self, marker, query_pair)


class Binder(ChromosomeSet):
    ''' A Binder is a special child class of ChromosomeSet which defines
    whether a genes is expressed as Left Maternal (inherit from Mothers
    left side and Fathers right, we denote this as 'LM') or Right Maternal
    (inherit from Mothers right side and Fathers left side, we denote this as
    'RM'). For Future references it might be useful to rename this class as
    Human_Binder as this only works class only works with Humans
    '''
    def __init__(self):
        '''(Binder) -> NoneType
        Initialize the Binder object when given a binder object
        '''
        # initialize the parent class with human chromosome with
        ChromosomeSet.__init__(self, human_chromosome, h_chromo_sex)
        # set sex to X until user sets the gender of the binder
        self._sex = 'X'

    def set_sex(self, sex):
        ''' (Binder, str) -> NoneType
        Sets the sex for of the Binder
        '''
        self._sex = sex

    def get_sex(self):
        ''' (Binder) -> NoneType
        Returns the sex of the binder
        '''
        return self._sex

    def set_by_pos(self, pair_num, pos, maternal):
        ''' (Binder, int, int, str) -> NoneType
        REQ: pair_num and pos >= 0
        REQ: maternal can only be maternal symbols by default these are
        specified to be 'LM' or 'RM'
        RAISE: ChromosomeNotFoundError
        Sets a gene pointed at by pair and position number to be a Left
        Maternal or Right Maternal.
        - Raises ChromosomeNotFoundError if ChromosomePair is not found
        '''
        # overwrite parent method for new requirements
        ChromosomeSet.set_by_pos(self, pair_num, pos, maternal)

    def set_by_marker(self, marker, maternal):
        ''' (Binder, str, str) -> NoneType
        REQ: maternal can only be maternal symbols by default these are
        specified to be 'LM' or 'RM'
        RAISE: MarkerNotFoundError
        Sets gene pointed at by Marker to be a Left Maternal or Right
        Maternal.
        - Raises MarkerNotFoundError if no marker was found in ChromosomeSet
        object
        '''
        ChromosomeSet.set_by_marker(self, marker, maternal)


class Male(Human):
    ''' A Male is a child of Human class
    '''
    def __init__(self, male_id):
        ''' (Human, int) -> NoneType
        Intialize a human object
        '''
        # initialize parent class and set it to have a set of chromosome
        # with human_chromosome number of chromosome
        Human.__init__(self, male_id)


class Female(Human):
    '''A female is a child of the Human class
    '''
    def __init__(self, female_id):
        ''' (Human, int) -> NoneType
        Intialize a female object
        '''
        # initialize parent class and set it to have a set of chromosome
        # with human_chromosome number of chromosome
        Human.__init__(self, female_id)

    def procreate(self, father, binder):
        ''' (Female, Male, Binder) -> obj
        RAISE: SexlessBabyError
        RAISE: InsufficientGeneError
        Return either Male or Female given a Female, Male and Binder. The
        offspring returned is determeind by the sex of the Binder. And each
        individual gene in the child is either Left Maternal (meaning inherit
        left nucloetide from Female and right from Male we denote this as
        'RM') or Right Maternal (inherit right nucleotide from Female and left
        from Male we denote this as 'LM'). Child will only inherit chromosome
        specified in the binder the other chromosomes pairs are assumed to be
        irrevelant.
        - Raise SexlessBabyError if Baby does not have a valid sex.
        - Raise InsufficientGeneError if either father or Mother
        is missing a gene that was specified in binder when procreating
        '''
        # define our global variables which are child born for our id of child
        # and left and right maternal symbols for binders and male and female
        # symbol to determine sex of the child
        global child_born, l_mat_sym, r_mat_sym, male_symbol, female_symbol
        # get babies sex using get sex method so we can determine babies sex
        sex = binder.get_sex()
        # if sex is not Male and Female
        if ((sex != male_symbol) and (sex != female_symbol)):
            # then our baby has an indeterminate sex and we cannot procreate
            raise SexlessBabyError('Binder sex is not ' + male_symbol +
                                   ' or ' + female_symbol)
        # otherwise define what our total chromosome is
        human_chromosome = self.get_total_pair()
        # if it is M
        if (sex == male_symbol):
            # create child to be Male with father id + mothers id + unique
            # identifier indicating how many child was born
            child = Male(father.get_id() + self.get_id() + str(child_born))
        # if it is F
        elif (sex == female_symbol):
            # create child to be Female with mothers id + father id + unique
            # identifier
            child = Female(self.get_id() + father.get_id() + str(child_born))
        # increment unique child born count to prevent a collision if couple
        # decides to procreate again
        child_born += 1
        # loop through all chromosome in binder since this a human binder
        # we loop from 0 to human_chromosome + 1 we add 1 to include
        # full range of the human chromosome
        for b_chromo_num in range(0, human_chromosome + 1):
            # get all position within current chromosome using get chromosome
            # and get all pos method
            b_chromo_pair = binder.get_chromosome(b_chromo_num).get_all_pos()
            # loop through each position
            for b_pos in b_chromo_pair:
                # try to access the positions in father and mother
                try:
                    mother_pair = self.get_by_pos(b_chromo_num, b_pos)
                    father_pair = father.get_by_pos(b_chromo_num, b_pos)
                    # if we get what is a the current index in binder and
                    # no error was raised get binder at current position
                    maternal = binder.get_by_pos(b_chromo_num, b_pos)
                    # get left and right from father
                    father_l = father_pair[0]
                    father_r = father_pair[1]
                    # get left and right from Mom
                    mother_l = mother_pair[0]
                    mother_r = mother_pair[1]
                    # check if it is Left Maternal
                    if (maternal is l_mat_sym):
                        # we get left from female right from male and set it to
                        # this position  in child
                        new_pair = mother_l + father_r
                        child.set_by_pos(b_chromo_num, b_pos, new_pair)
                    # if it is Right Maternal
                    elif (maternal is r_mat_sym):
                        # we get right from female left from male and set it to
                        # this position in child
                        new_pair = father_l + mother_r
                        child.set_by_pos(b_chromo_num, b_pos, new_pair)
                # if we do not find a gene in father or mother except the
                # Not Found Error
                except PairNotFoundError:
                    # raise the InsuffientGeneError because we must have
                    # gene in both parent otherwise we can't have offspring
                    raise InsuffientGeneError('Father or Mother' +
                                              ' have insufficient genes at' +
                                              ' pair: ' + str(b_chromo_num) +
                                              ' position: ' + str(b_pos))
        # return the child with all the genetic modification
        return child
