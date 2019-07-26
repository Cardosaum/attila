import collections
import pprint

cronograma_semanal = collections.defaultdict(dict)

for n in range(1, 18):
    cronograma_semanal[f'semana {n:02d}']['Matéria'] = ''
    cronograma_semanal[f'semana {n:02d}']['Teoria'] = ''
    cronograma_semanal[f'semana {n:02d}']['Exercícios'] = ''



def update_materia(week, content):
    cronograma_semanal[f'semana {week}']['Matéria'] = content

def update_teoria(week, content):
    cronograma_semanal[f'semana {week}']['Teoria'] = 'Tópicos: ' + content

def update_exercicios(week, content):
    cronograma_semanal[f'semana {week}']['Exercícios'] = content

update_materia('01', 'Revisão de Química Geral')
update_teoria('01', '1.1, 1.2, 1.3, 1.4, 1.5')

pprint.pprint(dict(cronograma_semanal))


