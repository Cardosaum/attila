from timeit import default_timer as timer

start = timer()


with open('/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq-TEST.txt', encoding='ISO-8859-1') as file:
    content = file.readlines()

print(content)

end = timer()

print(f'Elapsed time: {end - start}')
