import random, math, copy

with open("inputul.txt") as f:
    lines=f.read().splitlines()

dim_pop=int(lines[0])
dom_def=[int(x) for x in lines[1].split()]
coef_f=[int(x) for x in lines[2].split()]
precizie=int(lines[3])
p_recomb=float(lines[4])
p_mut=float(lines[5])
nr_epoch=int(lines[6])
ok=1
l_chrom=math.ceil(math.log((dom_def[1] - dom_def[0])*(math.pow(10,precizie)), 2))
g = open("evo.txt", "w")

#print(l_chrom)
def mutation(epoch):
	if epoch==1:
		g.write(f'Probabilitatea de mutatie pentru fiecare gena {p_mut}\nAu fost modificati cromozomii:\n')
	for i in range(dim_pop):
		u=random.random()
		if u<p_mut:
			if epoch==1:
				g.write(f'{i+1}\n')
			idx=random.randrange(0,l_chrom)
			pop_next[i][idx]=1-pop_next[i][idx]

def crossover(cross_idx, epoch):
	changed_ch1=[]
	changed_ch2=[]
	while len(cross_idx)>1:
		sample=random.sample(cross_idx,2)
		cross_idx.remove(sample[0])
		cross_idx.remove(sample[1])

		break_point=random.randrange(0,l_chrom+1)
		if epoch==1:
			g.write(f'Recombinare dintre cromozomul {sample[0]+1} cu cromozomul {sample[1]+1}:\n')
			g.write("".join([str(digit) for digit in pop_next[sample[0]]])+" "+"".join([str(digit) for digit in pop_next[sample[1]]]) +f' punct {break_point}\n')

		changed_ch1=pop_next[sample[0]][:break_point]+pop_next[sample[1]][break_point:]
		changed_ch2=pop_next[sample[1]][:break_point]+pop_next[sample[0]][break_point:]

		if epoch==1:
			g.write("Rezultat "+"".join([str(digit) for digit in changed_ch1])+" "+"".join([str(digit) for digit in changed_ch2])+"\n")
		pop_next[sample[0]]=changed_ch1
		pop_next[sample[1]]=changed_ch2

def cross_pop(pop, epoch):
	temp=[]
	if epoch==1:
		g.write(f'Probabilitatea de incrucisare {p_recomb}\n')
	for i in range(dim_pop):
		pp=random.random()
		afis=""
		afis+=f'{i+1}: '+"".join([str(digit) for digit in pop[i]])+f' u={pp}'
		if pp<p_recomb:
			temp.append(i)
			afis+=f'<{p_recomb} participa'
		if epoch==1:
			g.write(afis+"\n")
	return temp

def bs_interval(u, v, st, dr): 
    last=0
    while st <= dr:
        mij = (st + dr) // 2
        if v[mij] <= u:
            last = mij
            st = mij+1
        elif v[mij] > u:
            dr = mij-1
    return last

def selected_pop(intervale, epoch, pop):
	select=[]
	new_pop=copy.deepcopy(pop)
	for i in range(dim_pop):
		u=random.random()
		c_index=bs_interval(u,intervale,0,dim_pop)
		if epoch==1:
			g.write(f'u= {u} selectam cromozomul {c_index+1}\n')
		select.append(new_pop[c_index])
	return select

def interv_prob(pp, epoch):
	prob=[]
	x=0
	for i in pp:
		prob.append(x)
		x+=i
	prob.append(1.0)
	if epoch==1:
		g.write("Intervale probabilitati selectie\n")
		for i in prob:
			g.write(f'{i} ')
		g.write("\n")
	return prob

def prob_selectie(pop, epoch):
	all_fit=[func(x_getter(chrom)) for chrom in pop]
	fit_total=sum(all_fit)
	prob=[]
	for fit in all_fit:
		prob.append(fit/fit_total)
	if epoch==1:
		g.write("Probabilitati selectie\n")
		for i in range(dim_pop):
			g.write(f'Cromozom {i+1} probabilitate {prob[i]}\n')
		g.write("\n")
	return prob

def func(x):
	return coef_f[0]*(x**2)+coef_f[1]*x+coef_f[2]

def x_getter(chrom):
	x="".join([str(digit) for digit in chrom])
	x=int(x,2)
	return round(((dom_def[1]-dom_def[0]) / (2**l_chrom - 1))*x + dom_def[0], precizie)
def str_pop(pop):
	for i in range(dim_pop):
		x=x_getter(pop[i])
		f=func(x)
		g.write(f'{i+1}: '+"".join([str(digit) for digit in pop[i]])+f' x={x} f= {f}\n')
	g.write("\n")

def gen_chrom(l_chrom):
	return [random.randint(0,1) for j in range(l_chrom)]

def gen_pop():
	temp=[gen_chrom(l_chrom) for i in range(dim_pop)]
	return temp

def elite():
	idx1=0
	maxim=func(x_getter(pop[0]))
	for i in range(1,dim_pop):
		fit=func(x_getter(pop[i]))
		if fit>maxim:
			maxim=fit
			idx1=i
	idx2=0
	minim=func(x_getter(pop_next[0]))
	for i in range(1,dim_pop):
		fit=func(x_getter(pop_next[i]))
		if fit<minim:
			minim=fit
			idx2=i
	if epoch==1:
		g.write(f'Cromozomul elite e {idx1+1} din populatia initiala si inlocuieste {idx2+1}\n')

	pop_next[idx2]=copy.deepcopy(pop[idx1])

def val_max_mediu():
	suma=func(x_getter(pop_next[0]))
	maxim=func(x_getter(pop_next[0]))
	for i in range(1,dim_pop):
		fit=func(x_getter(pop_next[i]))
		suma+=fit
		if fit>maxim:
			maxim=fit
	suma=suma/dim_pop
	return maxim,suma

pop=gen_pop() #lista cu cromozomii din 1,0
#print(pop)
for epoch in range(nr_epoch):
	if epoch==1:
		g.write("Populatia initiala\n")
		str_pop(pop) #functie care afiseaza lista de chrom + x + f

	pp=prob_selectie(pop, epoch) #generare probabilitati pt fiecare si afis prob pt fiecare
									# cu formula f/sum(f)

	intervale_prob=interv_prob(pp, epoch) #returneaza lista cu intervalele suma probabilitati etc
	pop_next=selected_pop(intervale_prob, epoch, pop) #alege dim_pop prob random si adaugam obiectul asociat cu intervalul precedent de prob

	if epoch==1:
		g.write("Dupa selectie:\n")
		str_pop(pop_next)

	cross_idx_pop=cross_pop(pop_next, epoch) #idx cu chrom care au u<pc si participa in cross

	crossover(cross_idx_pop, epoch) #crossover pe o poz random cu elemente random din crossidxpop

	if epoch==1:
		g.write("Dupa recombinare:\n")
		str_pop(pop_next)

	mutation(epoch) #daca u<pm, pe acel chrom inversezi 1/0 pe un idx random
	if epoch==1:
		g.write("Dupa mutatie:\n")
		str_pop(pop_next)

	elite() #scoatem chrom cu f==min(f) din noua pop si adaugam f==max(f) din vechea pop
	if epoch==1:
		g.write("Dupa adaugare elite:\n")
		str_pop(pop_next)

	if epoch==1:
		g.write("Evolutia maximului\n")
	
	val_max, val_mediu=val_max_mediu() #returneaza f maxim din noua pop si media dintre toate f-urile
	if epoch>=1:
		g.write(f'f_max= {val_max} f_mediu= {val_mediu}\n')
	pop=copy.deepcopy(pop_next) #urmatoarea generatie
