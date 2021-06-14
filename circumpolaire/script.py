import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os
from dfply import *


class affichageDonnee:
	def __init__(self,fname):
		self.energieCinetique=list()
		self.fname=fname
		self.dict={'Potential_energy' : [],'Kinetic' : [],'Enstrophy' : [],'Transfer_potential_kinetic' : [],'PE_input_relax' : [],'PE_dissipation' : [],'KE_coriolis_nonlinear_terms' : [],'Transfer_kinetic_potential' : [],'KE_from_wind' : [],'KE_dissipation_bottom_friction' : [],'KE_dissipation_damping' : [],'KE_dissipation_lateral_diffusion' : []}
		self.dissipation=list()
		self.index=list()
		self.nbrLigne=0

	def RecupDonnee(self):
		self.AllDonnee=list()
		with open(self.fname,'r') as fichierdg:
			for l in fichierdg:
				if "#" in l:
					continue
				data=l.split()
				donnee=list()
				for i in range(len(data)):
					donnee.append(data[i])
					donnee[i]=float(donnee[i])
				self.AllDonnee.append(donnee)
				self.index.append(float(donnee[0]))
				self.energieCinetique.append(donnee[9])
				self.dissipation.append(donnee[12])
				self.dict['Potential_energy'].append(donnee[1])
				self.dict['Kinetic'].append(donnee[2])
				self.dict['Enstrophy'].append(donnee[3])
				self.dict['Transfer_potential_kinetic'].append(donnee[4])
				self.dict['PE_input_relax'].append(donnee[5])
				self.dict['PE_dissipation'].append(donnee[6])
				self.dict['KE_coriolis_nonlinear_terms'].append(donnee[7])
				self.dict['Transfer_kinetic_potential'].append(donnee[8])
				self.dict['KE_from_wind'].append(donnee[9])
				self.dict['KE_dissipation_bottom_friction'].append(donnee[10])
				self.dict['KE_dissipation_damping'].append(donnee[11])
				self.dict['KE_dissipation_lateral_diffusion'].append(donnee[12])
				self.nbrLigne+=1
		

	def graph(self):
		plt.plot(self.index,self.df.KE_from_wind,color='blue',label = 'Energie cinétique du vent donnée au fluide')
		plt.plot(self.index,self.df.KE_dissipation_lateral_diffusion,color='red', label ='Dissipation latérale du fluide')
		plt.plot(self.index,self.df.KE_dissipation_damping,color='black',label ='KE_dissipation_damping')
		plt.plot(self.index, self.df.KE_dissipation_bottom_friction,color='pink',label='KE_dissipation_bottom_friction')
		# plt.plot(self.index,self.df.Potential_energy,label ='Potential_energy')
		# plt.plot(self.index,self.df.Kinetic,label ='Kinetic')
		plt.plot(self.index,self.df.sumDis,label ='Somme dissipation')
		plt.xlabel('heures')
		plt.legend()
		plt.show()
		plt.savefig('graph_testvent_50000')

	def createDataFrame(self):
		
		self.df=pd.DataFrame(self.dict)
		self.df.index=self.df.index*1000/3600
		self.df.KE_dissipation_bottom_friction=-self.df.KE_dissipation_bottom_friction
		self.df.KE_dissipation_lateral_diffusion=-self.df.KE_dissipation_lateral_diffusion
		self.df.Kinetic=self.df.Kinetic/30
		

	def affichageDonnee(self):
		print(self.df.head())

	def sommeDissipation(self):
		self.df['sumDis']= self.df['KE_dissipation_bottom_friction']+self.df['KE_dissipation_lateral_diffusion']

	def createExcel(self):
		self.df.to_excel('Value.xlsx',index=False,sheet_name="xlambda_0.1")





if __name__ == '__main__':
	plt.close()
	test4=affichageDonnee('test4.dg')
	test4.RecupDonnee()
	test4.createDataFrame()
	test4.sommeDissipation()
	test4.graph()
	# test4.affichageDonnee()
	# test4.createExcel()
