import random
def Sample(set_var):
    return random.choice(tuple(set_var))


class IRule:
    INT32_MAX =  2^63 - 1

    def __init__(self, name: str, level: int, minimum_picks: int, maximum_picks: int):
        self.level = level

        self.name = name
        self.genomes = []

        self.variables = set()
        self.minimum_picks = minimum_picks
        self.maximum_picks = maximum_picks
        self.rolling_id = 0

    def AddGenome(self, genome):
        self.genomes.append(genome)

    def AddVariable(self, variable):
        self.variables.add(variable)

    def GetRollingId(self):
        self.rolling_id += 1
        return self.rolling_id


    @staticmethod
    def SplitIntoVarAndRule(string: str):
        variable_str = ""
        rule_str = ""

        tokens = string.split('(')

        variable_str = tokens[0]

        if len(tokens) == 2:
            if tokens[1].endswith(')'):
                rule_str = tokens[1][:-1]
            else:
                print("Expected ')'. Abort.")
                exit()

        return variable_str, rule_str

    @staticmethod
    def RuleFromString(rule_name: str, string: str):
        i_gl = 0
        for i in range(1, len(string)):
            if not string[:i].isdigit():
                break
            i_gl += 1

        if i > 0:
            numeric = int(string[:i])
            if i != len(string):
                if string[-1] == '+':
                    return IRule(rule_name , 0, numeric, IRule.INT32_MAX)
        return None

    def MustPick(self, remaining_picks):
        if len(self.variables) == self.maximum_picks:
            return Sample(self.picked)
        if len(self.variables) >= self.minimum_picks:
            return False
        #if len(self.variables)+remaining_picks == self.minimum_picks:
        #    return False


    def PickRandom(self):
        return random.choice(tuple(self.variables))

    def AddPick(self, pick: str):
        self.picked.add(pick)

    def ToString(self):
        return "{}[{},{}]".format(self.name, self.minimum_picks, self.maximum_picks)


class IVariable:
    def __init__(self, name: str):
        self.name = name
        self.value = None
        # if len(self.name) < 3: self.name = "XXX"
        self.genomes = []
        self.taxonomic_id = None
        self.rules = []
        self.level = -1

    def AddGenome(self, genome):
        self.genomes.append(genome)

    def SetValue(self, value: str):
        self.value = value

    def GetUnsetGenomes(self):
        return set([genome for genome in self.genomes if not genome.genome_var.GetValue()])

    def GetValue(self):
        return self.value

    def ToString(self):
        if self.GetValue(): return "{}({})".format(self.name, self.value)
        return self.name

    def ToString2(self):
        #if self.GetValue(): return "{}({}) level: ".format(self.name, self.value, self.level)
        #return self.name

        return "{}\t{}\t{}".format(self.name, self.level, self.value)

class IGenome:
    def __init__(self, name: str, lineage: str):
        size = len(lineage.split(';'))
        self.genome_var = None
        self.variables = [None]*size
        self.rules = [None]*size
        self.lineage = lineage
        self.name = name
        self.level = -1
        #self.Load(lineage)

    def SetRule(self, rule, level: int):
        self.rules[level] = rule

    def TaxonomyStr(self):
        names = [""] * len(self.variables)
        for i in range(len(self.variables)):
            if self.variables[i]: names[i] = self.variables[i].name
            elif self.rules[i]: names[i] = self.rules[i].name
        return ';'.join(names)

    def SetGenomeVar(self, genome_var):
        self.genome_var = genome_var

    def __eq__(self, other):
        for i in range(len(self.rules)):
            if self.rules[i] != other.rules[i]:
                return False
            if self.variables[i] != other.variables[i]:
                return False
        return True

    def __hash__(self):
        return hash((tuple(self.rules), tuple(self.variables)))

    def HasRule(self, rule: IRule):
        return rule in self.rules

    def SetVariable(self, variable, level: int):
        self.variables[level] = variable

    def ToString(self):
        # placeholder = '___'
        rules_str = ""
        variables_str = ""

        for i in range(len(self.rules)):
            rule = self.rules[i]
            variable = self.variables[i]

            if i > 0: 
                rules_str += '  '
                variables_str += '  '

            placeholder_len = 3
            # if rule: placeholder_len = len(rule.ToString())
            # if variable: placeholder_len = max(placeholder_len, len(variable.ToString()))
            # if rule and variable: 
            #     placeholder_len = max(len(rule.ToString()), len(variable.ToString()))
            #     print(placeholder_len)

            placeholder_len = max(len(rule.ToString()) if rule else placeholder_len, len(variable.ToString()) if variable else placeholder_len)
            placeholder = '_' * placeholder_len

            vstr = placeholder
            rstr = placeholder
            if rule: 
                rstr = rule.ToString()
                if len(rstr) < placeholder_len:
                    rstr += ' ' * (placeholder_len - len(rstr))
            if variable:
                vstr = variable.ToString()
                if len(vstr) < placeholder_len:
                    vstr += ' ' * (placeholder_len - len(vstr))

            
            variables_str += "{}".format(vstr)
            rules_str += "{}".format(rstr)

        return "{} ({}):\n\tvariables: {}\n\trules:     {}".format(self.name, self.genome_var.GetValue(), variables_str, rules_str)


class INode:
    def __init__(self, variable, rule):
        self.variable = variable
        self.rule = rule
        self.parent = None
        self.children = set()

    def SetParent(self, parent):
        self.parent = parent

    def AddChild(self, child):
        self.children.add(child)

class Group:
    def __init__(self):
        self.members = []
        self.select = []

    def AddGenome(self, genome):
        self.members.append(genome)

    def SetSelect(self, select_str: str):
        if select_str.isdigit():
            self.select.append(int(select_str))

    def Select(self):
        select = random.choice(tuple(self.select))
        indices = random.sample(range(0, len(self.members)), select)
        return [self.members[i] for i in indices]

class Instructions:
    VARIABLE_HEADER = "#Variables"
    GENOME_COLUMN = 0
    LINEAGE_COLUMN = 1
    GROUP_COLUMN = 2
    SELECT_COLUMN = 3

    GENOMES_HEADER = "#Genomes"
    VARIABLE_COLUMN = 0
    RELATIVE_ABUNDANCE_COLUMN = 1

    def __init__(self, path: str):
        self.genomes = set()
        self.variables = dict()
        self.rules = dict()
        self.groups = dict()
        self.Load(path)

        self.rolling_gid = 0


    @staticmethod
    def IsComment(line : str):
        return line.startswith('#')

    @staticmethod
    def IsVariableHeader(line: str):
        return line == Instructions.VARIABLE_HEADER

    @staticmethod
    def IsGenomesHeader(line: str):
        return line == Instructions.GENOMES_HEADER

    def GetRollingGid(self):
        self.rolling_gid += 1
        return self.rolling_gid

    def GetVariablesOfLevel(self, level: int):
        return [var for var in self.variables.values() if var.level == level]

    def GetRulesOfLevel(self, level: int):
        return [var for var in self.rules.values() if var.level == level]

    def GetValues(self):
        return [var.GetValue() for var in self.variables.values() if var.GetValue()]

    def AddGroup(self, group_id: int) -> Group:
        if group_id not in self.groups:
            group = Group()
            self.groups[group_id] = group
        return self.groups[group_id]

    def AddRule(self, variable_str: str, rule_str: str) -> IRule:
        # print("Rule:     {}".format(rule.ToString()))
        if variable_str not in self.rules:
            rule = IRule.RuleFromString(variable_str, rule_str)
            self.rules[variable_str] = rule
        return self.rules[variable_str]

    def AddVariable(self, variable_str: str) -> IVariable:
        if variable_str not in self.variables:
            self.variables[variable_str] = IVariable(variable_str)
        return self.variables[variable_str]

    def GetEligibleChoices(self, taxonomy, genome_subset, level: int, name: str):
        genome_subset = list(genome_subset)
        genome_subset.sort(key=lambda g: g.TaxonomyStr())

        non_eligible = []

        print("GetEligibleChoices({}, {}, {})".format(level, name, genome_subset))

        ranks = len(genome_subset[0].variables)

        rank_requirements = dict()

        for genome in genome_subset:
            next_set = level + 1
            # if next_set == 

            if next_set < ranks:
                while not genome.variables[next_set]: next_set += 1

            if next_set not in rank_requirements: rank_requirements[next_set] = set()

            rname = genome.name if next_set == ranks else genome.variables[next_set].name
            rset = rank_requirements[next_set]
            rset.add(rname)
            
        # maybe provide as variable

        nodes = taxonomy.LeafsFromNode(name) if level == ranks else taxonomy.NodesFromNode(name, level+1)
        og_len = len(nodes)
        nodes = set(nodes).difference(set(node.GetValue() for node in self.GetVariablesOfLevel(level) if node.GetValue()))

        print("=>>>> {} -> {}".format(og_len, len(nodes)))

        print("These are the choices: From {} nodes: {}".format(name, len(nodes)))
        
        print("Nodes: {}".format(nodes))
        for rank,value in rank_requirements.items():
            min_requirement = len(value)
            print("rank: {} needs {} options".format(rank, value))

            # Check all nodes
            for node in nodes:
                # Iterate all potential choices for given variable (say, phylum)
                # Only for this node.
                nset = taxonomy.LeafsFromNode(node) if rank == ranks else taxonomy.NodesFromNode(node, rank+1)
                nset_len_og = len(nset)
                nset = set(nset).difference(set([var.GetValue() for var in self.GetVariablesOfLevel(rank)]))
                print("Charlie {}: {} -> {}".format(node, nset_len_og, len(nset)))
                print("Charlie {} -> {} (req: {}) Leaves? {}".format(node, len(nset), min_requirement, (rank == ranks)))
                
                if len(nset) < min_requirement:
                    non_eligible.append(node)
                    print("DROP: {}  non-eligible {}".format(node, non_eligible))
                else:
                    if rank == ranks: 
                        # print("Mach weiter weil species.")
                        continue

                    # of eligible nodes
                    for var in value:
                        # check each variable under the the current one
                        enum = 0
                        tested = 0
                        # for newnode in nset:
                            # remove the already chosen ones
                        new_genome_sub = self.variables[var].GetUnsetGenomes()
                        if len(new_genome_sub) == 0: continue
                        tested += 1
                        print("\n\n\nRECURSION: ################ var {}, target rank {}, old lca {}, len new genome_sub {}".format(var, rank, node, len(new_genome_sub)))
                        if len(self.GetEligibleChoices(taxonomy, new_genome_sub, rank, node)) > 0:
                            enum += 1
                        # if tested > 0 and enum == 0:
                           # non_eligible.append(node)

        print("Old set: {}".format(len(nodes)))
        print("Drop nodes: {}".format(len(non_eligible)))
        print("New set: {}".format(len(set(nodes).difference(set(non_eligible)))))
            # nodes = taxonomy.NodesFromNode("g__Phocaeicola", 7)
        
        eligible = set(nodes).difference(set(non_eligible))

        print("Remaining: {}".format(eligible))

        return eligible




    def ProcessGenome(self, line: str):
        tokens = line.split('\t')
        genome = tokens[Instructions.GENOME_COLUMN]
        lineage = tokens[Instructions.LINEAGE_COLUMN]
        igenome = IGenome(genome, lineage)


        genome_var = self.AddVariable(genome)
        genome_var.level = len(lineage.split(';'))
        igenome.SetGenomeVar(genome_var)


        # Process group column
        group_id = None
        if Instructions.GROUP_COLUMN < len(tokens):
            group_id = int(tokens[Instructions.GROUP_COLUMN])
        else:
            group_id = self.GetRollingGid()
        group = self.AddGroup(group_id)
        group.AddGenome(igenome)

        # Process group select column
        if Instructions.SELECT_COLUMN < len(tokens):
            group.SetSelect(tokens[Instructions.SELECT_COLUMN])

        self.genomes.add(igenome)

        lineage_tokens = lineage.split(';')

        for i in range(len(lineage_tokens)):
            lin = lineage_tokens[i]
            lin = lin.split('__')[-1] if "__" in lin else lin
            variable_str, rule_str = IRule.SplitIntoVarAndRule(lin)

            if variable_str or rule_str:
                print("Variable: {}\t\t{}".format(variable_str, rule_str))

            if rule_str:
                rule = self.AddRule(variable_str, rule_str)
                igenome.SetRule(rule, i)
                rule.AddGenome(igenome)
                rule.level = i

            elif variable_str:
                variable = self.AddVariable(variable_str)
                igenome.SetVariable(variable, i)
                variable.AddGenome(igenome)
                variable.level = i


    def Load(self, path: str):
        genome_section = False
        variable_section = False

        with open(path, 'r') as file:
            expect_header=True
            for line in file:
                line = line.rstrip()
                # print(line)

                if Instructions.IsGenomesHeader(line): 
                    print("genome section")
                    genome_section = True
                    variable_section = False
                    expect_header = True
                    continue
                elif Instructions.IsVariableHeader(line): 
                    print("variable section")
                    variable_section = True
                    genome_section = False
                    expect_header = True
                    continue

                elif Instructions.IsComment(line): continue

                tokens = line.split('\t')
                if genome_section:
                    if expect_header:
                        expect_header = False
                        continue

                    self.ProcessGenome(line)

                elif variable_section:
                    if expect_header:
                        expect_header = False
                        continue

                    print(line)


class GenomeSelector:
    @staticmethod
    def Select(instructions: Instructions, taxonomy, genomes):
        for genome in instructions.genomes:
            print(genome.ToString())
            GenomeSelector.SelectGenome(instructions, genome, taxonomy, genomes)

        for genome in instructions.genomes:
            print(genome.ToString())

    @staticmethod
    def SelectGenome(instructions, igenome: IGenome, taxonomy, genomes):
        print(igenome.ToString())

        # # Find first set from right
        # for i in range(len(igenome.variables)-1, -1, -1):
        #     if not igenome.variables[i] and not igenome.rules[i]:
        #         continue

        #     variable = igenome.variables[i]
        #     rule = igenome.rules[i]

        #     if variable and not variable.GetValue():
        #         print("Not set yet")

        genome_subset = set(genomes)

        last_set = taxonomy.RootName()

        # First unset from left
        for i in range(len(igenome.variables)):
            level = i+1
            print("__{}\t{}:".format(i, level))

            variable = igenome.variables[i]
            rule = igenome.rules[i]

            if variable:
                if variable.GetValue():
                    last_set = variable.GetValue()
                    continue

                print("Unset Variable: (last set {})".format(last_set))
                nodes = set(taxonomy.NodesFromNode(last_set, level))
                nodes = nodes.difference(set([var.GetValue() for var in instructions.GetVariablesOfLevel(i)]))
                print("pick from: {}".format(len(nodes)))

                nodes = instructions.GetEligibleChoices(taxonomy, variable.genomes, i, last_set)

                # input("GetEligibleChoices:")


                sampled = Sample(nodes)
                print("Pick: {}".format(sampled))
                print("Assume everything is fine......")

                variable.SetValue(sampled)
                last_set = variable.GetValue()

                print(igenome.ToString())
                #i = input("Picked: {}".format(sampled))

            elif igenome.rules[i]:
                print("Unset rule")


                remaining_picks = set(rule.genomes)
                print(remaining_picks)

                pick = rule.MustPick(len(remaining_picks))

                if pick:
                    random_var = pick.PickRandom()

                else:
                    nodes = taxonomy.NodesFromNode(last_set, level)
                    nodes = set(nodes).difference(instructions.GetValues()).union(set(var.GetValue() for var in rule.variables))

                    print("pick random freely")
                    print("pick from: {}".format(nodes))

                    sampled = Sample(nodes)
                    print("Pick: {}".format(sampled))
                    newvar_name = "{}_{}".format(rule.name, rule.GetRollingId())

                    newvar = None
                    if [var for var in instructions.variables.values() if var.GetValue() == sampled]:
                        newvar_name = [var for var in instructions.variables.values() if var.GetValue() == sampled][0].name
                        newvar = instructions.variables[newvar_name]
                        print("Exists: {} -> {}".format(sampled, newvar_name))
                    else:
                        newvar = IVariable(newvar_name)
                        newvar.SetValue(sampled)
                        instructions.variables[newvar_name] = newvar
                        newvar.level = i
                        rule.AddVariable(newvar)
                    
                    igenome.variables[i] = newvar
                    last_set = newvar.GetValue()

                    for j in range(i+1, len(igenome.variables)):
                        if igenome.variables[j]:
                            for genome in igenome.variables[j].genomes:
                                print("Update {}".format(genome.ToString()))
                                genome.variables[i] = newvar
                            break

                print(igenome.ToString())
                #i = input("input...")
            else:
                print("Nothing")

        genome_set = set(taxonomy.LeafsFromNode(last_set))
        print("Old set: {}".format(genome_set))
        genome_set = genome_set.difference(set([var.GetValue() for var in instructions.GetVariablesOfLevel(len(igenome.variables))]))
        print("New set: {}".format(genome_set))
        
        print("{} Existing leafs: {}".format(igenome.genome_var.name, ','.join(map(lambda x: str(x.GetValue()), instructions.GetVariablesOfLevel(len(igenome.variables))))))
        # for var in instructions.variables.values():
        #     print(var.ToString2())

        #input()
        genome_var = instructions.variables[igenome.name]
        choice = Sample(genome_set)
        genome_var.SetValue(choice)

        print("__________________ LAST SET {}".format(last_set))
        print("select from {}".format(genome_set))
        print("SELECT:  {}".format(choice))