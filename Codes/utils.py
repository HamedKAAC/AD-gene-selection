import pandas as pd
import numpy as  np


def load_network(path_to_nodes, path_to_edges, path_to_platform):
  net_nodes = pd.read_csv(path_to_nodes)
  net_edges = pd.read_csv(path_to_edges)

  node_ensp2name = {}
  node_name2ensp = {}
  for i in range(len(net_nodes)):
      node_ensp2name[net_nodes.iloc[i]['stringdb::database identifier']] = net_nodes.iloc[i]['display name']
      node_name2ensp[net_nodes.iloc[i]['display name']] = net_nodes.iloc[i]['stringdb::database identifier']

  coexp_net_edges = []
  for i in range(len(net_edges)):
      if net_edges.iloc[i]['interaction'] == 'pp' and net_edges.iloc[i]['stringdb::coexpression'] != 'nan':
          if float(net_edges.iloc[i]['stringdb::coexpression'])>.5:
            coexp_net_edges.append(net_edges.iloc[i]['name'])

  platform = open(path_to_platform, 'r').readlines()
  common_genes = []
  for line in platform[17:]:
    line_content = line.split('\t')
    if line_content[10] != '':
        gene_symbols = line_content[10].split('///')
        gene_symbols_available_in_network = []
        for gene_symbol in gene_symbols:
            if gene_symbol in node_name2ensp.keys():
                gene_symbols_available_in_network.append(gene_symbol)
        if len(gene_symbols_available_in_network)>=1:
            common_genes.append(gene_symbols_available_in_network[0])

  common_genes = list(set(common_genes))
  gene_adj_index = {}
  for i in range(len(common_genes)):
      gene_adj_index[common_genes[i]] = i
  adj = np.zeros((len(common_genes), len(common_genes)))
  for edge in coexp_net_edges:
      edge = edge.split('(pp)')
      from_ = edge[0].strip()
      to = edge[1].strip()
      if node_ensp2name[from_] in gene_adj_index.keys():
          from_ = gene_adj_index[node_ensp2name[from_]]
      else:
          continue
      if node_ensp2name[to] in gene_adj_index.keys():
          to = gene_adj_index[node_ensp2name[to]]
      else:
          continue
      adj[from_, to] = 1
      adj[to, from_] = 1

  return adj, gene_adj_index, common_genes


def Analysis(name_of_group, DEG_Data, degs_threshold, adj, gene_adj_index, common_genes):

  Group_DEG_Data = DEG_Data[:degs_threshold] 
  del Group_DEG_Data['Gene.title']
  Group_DEG_Data['evidence_score'] = [float('nan') for i in range(len(Group_DEG_Data))]


  def is_DEG(g):
    for deg in Group_DEG_Data['Gene.symbol']:
      if str(deg) != 'nan':
        if g in deg:
          return True
    return False

  for i in range(int(degs_threshold/2)):

    gene_symbols = Group_DEG_Data['Gene.symbol'][i]
    if str(gene_symbols) != 'nan':
      gene_symbols = gene_symbols.split('///')
      used_symbols_in_network = []
      for gene in gene_symbols:
        if gene in common_genes:
          used_symbols_in_network.append(gene)
      if len(used_symbols_in_network)>=1:
        focused_gene = used_symbols_in_network[0]

        focused_gene_index_in_adj = gene_adj_index[focused_gene]
        neighbors = np.where(adj[focused_gene_index_in_adj, :] == 1)[0]
        if len(neighbors) == 0:
 
          Group_DEG_Data.loc[i, 'evidence_score'] = 0.5
        else:
          number_of_neighbors_presented_in_DEG = 0
          for n in neighbors:
            if is_DEG(common_genes[n]):
              number_of_neighbors_presented_in_DEG += 1

          Group_DEG_Data.loc[i, 'evidence_score'] = (number_of_neighbors_presented_in_DEG/len(neighbors))
      else:
          Group_DEG_Data.loc[i, 'evidence_score'] = 0.5

  return Group_DEG_Data

def filter_evidence_score(genes_with_evidence_score):
  selected_genes = []
  for i in range(len(genes_with_evidence_score)):
    if genes_with_evidence_score.iloc[i]['evidence_score'] >= 0.5:
      selected_genes.append(i)

  return genes_with_evidence_score.loc[selected_genes]


def selected_genes(group_data, evidence_score_threshold, adj, gene_adj_index, common_genes):
  group_selected_genes = []
  all_genes = set()
  for gd in group_data:
    gdd = Analysis(gd[0], gd[1], evidence_score_threshold, adj, gene_adj_index, common_genes)
    filteredgdd = filter_evidence_score(gdd)
    group_selected_genes.append(
        list(set(filteredgdd['Gene.symbol'].values.tolist()))
    )
    all_genes = all_genes.union(set(filteredgdd['Gene.symbol'].values.tolist()))



  all_genes = list(all_genes)
  genes_count = []
  for gene in all_genes:
    gene_counter = 0
    for gsgs in group_selected_genes:
      if gene in gsgs:
        gene_counter += 1
    genes_count.append(gene_counter)
  common_in_at_least_three = []
  for i in range(len(genes_count)):
    if genes_count[i]>= 3:
      common_in_at_least_three.append(
          all_genes[i]
      )


  last_group_selected_genes = []
  for gsgs in group_selected_genes:
    temp = []
    for gene in gsgs:
      if gene in common_in_at_least_three:
        temp.append(
            gene
        )
    last_group_selected_genes.append(temp)

  return common_in_at_least_three, last_group_selected_genes


def load_data(path_to_data, path_to_labels, path_to_platform, additionalMappings = {}):
  mappings = {}
  with open(path_to_platform, 'r') as map_file:
    lines = map_file.readlines()
    for line in lines[17:]:
      line_content = line.split('\t')
      if not line_content[10] == '':
        mappings[line_content[0].replace(' ', '')] = line_content[10].replace(' ', '')

  mappings = {**mappings, **additionalMappings}


  Data = pd.read_csv(path_to_data, index_col = 0, sep = '\t').T[list(mappings.keys())]
  Data = Data.rename(columns=mappings)
  unique_columns = list(set(Data.columns.tolist()))
  newData = pd.DataFrame()
  for col in unique_columns:
    temp = Data[[col]]
    newData[col] = temp.sum(axis=1)/temp.shape[1]
  newData = np.log2(newData)
  for col in newData.columns:
      newData[col] = (newData[col] - newData[col].mean())/newData[col].std(ddof=0)

  Labels = pd.read_csv(path_to_labels, usecols = ['GEO Accession', 'Disease State'] , index_col=0)
  y = []
  for i in range(161):
    if 'normal' in Labels.iloc[i]['Disease State']:
      y.append(1)
    else:
      y.append(0)

  return newData, y, Labels