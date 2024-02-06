from astropy.io import fits
import numpy as np
import math
import matplotlib.pyplot as plt

chandra = fits.open('chandra.fits') # Enter your fits file here. 
chandra.info()
data = chandra[1].data

print("--------------------------------------------PRIMARY HDU HEADER--------------------------------------------")
print(chandra[0].header)
print("--------------------------------------------PRIMARY HDU DATA--------------------------------------------")
print(chandra[0].data)
print("--------------------------------------------BINTABLE HDU 1 (EVENTS) DATA-------------------------------------")

print("--------------------------------------------BINTABLE HDU 2 DATA--------------------------------------------")
data2 = chandra[2].data


print("--------------------------------------------BINTABLE HDU 1 COLUMNS--------------------------------------------")
print(chandra[1].columns)
print("--------------------------------------------BINTABLE HDU 1 ENERGIES--------------------------------------------")
chandraData1 = chandra[1].data
print(chandraData1['energy'])
print(chandraData1['energy'][0:5])

# energy unit is eV, keV= eV/1000
print("-------------------------------------------FILTERED ENERGIES --------------------------------------------")
chandraData1['energy']/=1000
pha_values = chandraData1['energy']
x_values = chandraData1['X']
y_values = chandraData1['Y']

# Use numpy.where to filter the data based on energy range criteria
energy_mask = np.where((pha_values >= 1.5) & (pha_values <= 6))
energies = pha_values[energy_mask]
x_positions = x_values[energy_mask]
y_positions = y_values[energy_mask]
result = np.concatenate((energies.reshape(-1,1),x_positions.reshape(-1,1),y_positions.reshape(-1,1)),axis =1)
print("Now results: ")
print(result)
result_copy=result.copy()


#---------------------------------Pixel Merging to find center----------------------------------------------
merged_results = {}

for subarray in result:
    energy, x, y = subarray
    
    if int(x)<=x<int(x)+0.5:
        x = int(x)
    else:
        x = int(x)+1
    
    if int(y)<=y<int(y)+0.5:
        y = int(y)
    else:
        y = int(y)+1

    key = (x, y)

    if key in merged_results:
        merged_results[key] += energy
    else:
        merged_results[key] = energy

merged_list = [(energy, x, y) for (x, y), energy in merged_results.items()]
merged_array = np.array(merged_list)

class point:
    def __init__(self, x, y, energy):
        self.x = x
        self.y = y
        self.energy = energy


max_energy_index = np.argmax(merged_array[:, 0])

center_energy = merged_array[max_energy_index, 0]
center_x = merged_array[max_energy_index, 1]
center_y = merged_array[max_energy_index, 2]

center_point = point(center_x, center_y, center_energy)

print("Center Point:")
print("X:", center_point.x)
print("Y:", center_point.y)
print("Energy:", center_point.energy)
#--------------------------------------------------------------------------------------------------------------

pixel_arcsecond=0.492
maxRadius = 350 #arcsec
minRadius = 20
interval= 5


#---------------------------------Convertion - Pixel to Arcsecond----------------------------------------------
def convert():
    for sublist in result: #now all pixels are in arcsecond in result
        sublist[1]=sublist[1]*pixel_arcsecond
        sublist[2]=sublist[2]*pixel_arcsecond

#--------------------------------------------------------------------------------------------------------------




#---------------------------------Calculates Photon Intensity Within Chip----------------------------------------

def calculate_area(rad1, rad2): #in arcsec
    return math.pi*(rad2**2-rad1**2)

radius_intensity={}
def photon_intensity(minRadius, maxRadius, center):
    convert()
    for i in range(minRadius, maxRadius+1, 5):
        radius_intensity[i]=0
    
    for point in result:
        x=point[1]
        y=point[2]
        distance = math.sqrt(((x - center.x*pixel_arcsecond) ** 2 + (y - center.y*pixel_arcsecond) ** 2))
        
        if distance<=maxRadius:
            if distance<minRadius:
                radius_intensity[minRadius]+=1
            elif distance%5!=0:
                radius_intensity[int(distance+(5-distance%5))]+=1
            else:
                radius_intensity[distance]+=1

    for radius in radius_intensity.keys():
        if radius==minRadius:
            radius_intensity[radius]=radius_intensity[radius]/calculate_area(0,radius)
        else:
            radius_intensity[radius]=radius_intensity[radius]/calculate_area(radius-5,radius)
#----------------------------------------------------------------------------------------------------------------





#---------------------------------Energy Intensity----------------------------------------------
photonIntensities = {}

#Create 0.1 wide intervals
for i in range(75, 301):
    photonIntensities[i/50] = 0

#Round energies to their 1/100 decimals, match them to corresponding energy levels within
for i in result_copy:
    if pixel_arcsecond*(((i[1]-center_point.x)**2+(i[2]-center_point.y)**2)**(0.5)) < maxRadius:
        energy = round(i[0]*50)/50
        for j in photonIntensities.keys():
            if j==energy:
                photonIntensities[j] += 1    

intensityX= []
intensityY = []
for i, j in photonIntensities.items():
    intensityX.append(i)
    intensityY.append(j)
#----------------------------------------------------------------------------------------------






#---------------------------------Calculates Photon Intensity Within and Out of Chip-----------------------------

coordinatesEnergy = {}
for i in range(minRadius, maxRadius+1, interval):
    coordinatesEnergy[i] = 0 # The key refers to the positioning of photons between the key's circle and previous circle
    


#Calculating photon count per arcsec^2 
for i in result_copy:
    radius = pixel_arcsecond*(((i[1]-center_point.x)**2+(i[2]-center_point.y)**2)**(0.5))
    if (radius <= maxRadius):
        for j in coordinatesEnergy.keys():
            if radius > j:
                pass
            else:
                coordinatesEnergy[j] += 1
                break
        

countPerArcsecSq =  {}
prevRadius = 0
intervalCounter = interval
edgeCircle = 258
prevPieMinusTri = 0
for i, j in coordinatesEnergy.items():
    currentRadius = i
    
    if (currentRadius > edgeCircle):
        
        degrees = ( 180.0 / math.pi ) * math.acos(edgeCircle/(edgeCircle+intervalCounter))
        areaPart1 = math.pi*((edgeCircle+intervalCounter)**2)*(degrees/90)
        areaPart2 = math.sin(2*math.acos(edgeCircle/(edgeCircle+intervalCounter)))*((edgeCircle+intervalCounter)**2)
        pieMinusTri = areaPart1 - areaPart2
        emptyArea = pieMinusTri - prevPieMinusTri
        intervalRegion = (math.pi*(currentRadius)**2)- (math.pi*(prevRadius)**2)-emptyArea
       
        countPerArcsecSq[i] = (j / intervalRegion)
        prevRadius=currentRadius
        intervalCounter += interval
        prevPieMinusTri = pieMinusTri
        
        
    else:
        area = (math.pi*currentRadius**2)-(math.pi*prevRadius**2)
        countPerArcsecSq[i] = (j / area)
        prevRadius = currentRadius
        
#----------------------------------------------------------------------------------------------------------------      


#This part prints the chart with radius values within and out of the chip

x = list(countPerArcsecSq.keys())
y = list(countPerArcsecSq.values())

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))

ax1.plot(intensityX, intensityY, color="darkgreen")
ax2.plot(x, y, color="gold")

ax1.set_title("Photon Count Within " + str(maxRadius) + " Arcsec Radius (1.5, 6 keV)")
ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel('Photon Count')


ax2.set_title('Photon Intensity' + " (" + str(1.5) + ", " + str(6) + " keV)")
ax2.set_xlabel('Distance (Arcsec)')
ax2.set_ylabel('Intensity (count/arcsec^2)')


plt.tight_layout()  # Improves spacing between elements
plt.show()



#This part prints the chart with radius values within the chip
"""
photon_intensity(minRadius, maxRadius, center_point)

x_values = list(radius_intensity.keys())
y_values = list(radius_intensity.values())

plt.plot(x_values, y_values, color='blue', linewidth=2, markersize=6)
plt.xlabel('Radius (arcsec^2)')
plt.ylabel('Intensity (photon count/arcsec^2)')
plt.title('Photon intensity vs Radius (1.5 - 6.0 keV)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xticks(range(0, max(x_values)+1, 25), fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()  # Improves spacing between elements
plt.show()
"""