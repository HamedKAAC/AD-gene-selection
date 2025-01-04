import random
import math
from sklearn.neighbors import KNeighborsClassifier
from tqdm import tqdm

def initialize_poplation(L):
  population = []
  n = math.floor(1.65*2**(0.21*L))
  if n%2 == 1:
    n = n + 1

  for i in range(n):
    chromosome = {
        'structure': [],
        'fitness': None
    }

    for j in range(L):
      if random.uniform(0, 1) > 0.5:
        chromosome['structure'].append(1)
      else:
        chromosome['structure'].append(0)

    population.append(chromosome)

  return population





def calculate_fitness(population, dataframe, labels, save):

  features = dataframe.columns.values.tolist()
  for chromosome in population:
    selected_features = []
    for i, status in enumerate(chromosome['structure']):
      if status == 1:
        selected_features.append(features[i])

    try:
      fitness = save[str(chromosome['structure'])]
      chromosome['fitness'] = fitness
    except:
      if len(selected_features) == 0:

        chromosome['fitness'] = 0
        save[str(chromosome['structure'])] = chromosome['fitness']
      else:

        X = dataframe[selected_features].values.tolist()

        y = labels

        c = 0
        for i in range(len(dataframe)):
          X_train, X_test = X[0:i]+X[i+1:], [X[i]]
          y_train, y_test = y[0:i]+y[i+1:], [y[i]]
          clf = KNeighborsClassifier(3)
          clf.fit(X_train, y_train)
          if clf.score(X_test, y_test) == 1:
            c+=1

        n = sum(chromosome['structure'])
        chromosome['fitness'] = 0.8*c+0.2/n
        save[str(chromosome['structure'])] = chromosome['fitness']


  return population


def selection(population):

    avg = sum([chromosome['fitness'] for chromosome in population])/len(population)

    proportions = [chromosome['fitness']/avg for chromosome in population]

    number_of_copies = []
    for prop in proportions:

      int_ = int(prop)

      float_ = prop - int_

      rand = random.uniform(0, 1)

      if float_ > rand:
        number_of_copies.append(int_ + 1)
      else:
        number_of_copies.append(int_)

    mating_pool = []
    for i in range(len(population)):
      for j in range(number_of_copies[i]):
            mating_pool.append(population[i])

    if len(mating_pool) < len(population):
      difference = len(population) - len(mating_pool)
      for i in range(difference):
        mating_pool.append(population[-1])

    return mating_pool


def crossover(population):
  new_population = []
  
  if len(population)%2==1:
    iter = len(population)-1
  else:
    iter = len(population)
  
  for i in range(0, iter, 2):
    if random.uniform(0, 1) < .9:
      p1 = population[i]['structure']
      p2 = population[i+1]['structure']
      rand = random.randint(0, len(p1))
      ch1 = p1[0:rand] + p2[rand:]
      ch2 = p2[0:rand] + p1[rand:]
      new_population.append(
            {'structure': ch1, 'fitness': None}
        )
      new_population.append(
            {'structure': ch2, 'fitness': None}
      )
    else:
      p1 = population[i]['structure']
      p2 = population[i+1]['structure']
      new_population.append(
            {'structure': p1, 'fitness': None}
        )
      new_population.append(
            {'structure': p2, 'fitness': None}
      )

  return new_population

def mutation(population):
  new_population = []
  for chromosome in population:
    chromosome_structure = chromosome['structure']
    if random.uniform(0, 1) < 0.05:
      rand = random.randint(0, len(chromosome_structure)-1)
      if chromosome_structure[rand] == 0: chromosome_structure[rand] =1
      else: chromosome_structure[rand] = 0

    new_population.append({'structure': chromosome_structure,'fitness': None})

  return new_population

def best(elite_chromosome, population):
  best_fitness = elite_chromosome['fitness']
  best_chromosome = elite_chromosome
  for chromosome in population:
    if chromosome['fitness'] >= best_fitness:
      best_fitness = chromosome['fitness']
      best_chromosome = chromosome

  return best_chromosome
def average(population):
  return sum([ch['fitness'] for ch in population])/len(population)

def GA(dataframe, labels, iteration):
    save = {}

    L = len(dataframe.columns.values.tolist())

    population = initialize_poplation(L)

    population = calculate_fitness(population, dataframe, labels, save)

    elitism_variable = best({'structure': None, 'fitness': float('-inf')}, population)

    for iter in range(iteration):
      print(iter, elitism_variable)
      population = selection(population)
      population = crossover(population)
      population = mutation(population)
      population[random.randint(0, len(population)-1)] = elitism_variable
      population[random.randint(0, len(population)-1)] = elitism_variable
      population = calculate_fitness(population, dataframe, labels, save)
      elitism_variable = best({'structure': None, 'fitness': float('-inf')}, population)

    return elitism_variable