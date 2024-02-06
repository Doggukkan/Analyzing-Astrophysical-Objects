from astropy.io import fits
import numpy as np
import math
import matplotlib.pyplot as plt

chandra = fits.open('swift.evt') # Enter your fits file here
chandra.info()
data = chandra[1].data
array = data.view(np.recarray)
print("--------------------------------------------PRIMARY HDU HEADER--------------------------------------------")
print(chandra[0].header)
print("--------------------------------------------PRIMARY HDU DATA--------------------------------------------")
print(chandra[0].data)
print("--------------------------------------------BINTABLE HDU 1 (EVENTS) DATA-------------------------------------")
print(array)
print("--------------------------------------------BINTABLE HDU 2 DATA--------------------------------------------")
data2 = chandra[2].data
array2 = data2.view(np.recarray)
print(array2)
print("--------------------------------------------BINTABLE HDU 1 COLUMNS--------------------------------------------")
print(chandra[1].columns)
print("--------------------------------------------BINTABLE HDU 1 ENERGIES--------------------------------------------")
chandraData1 = chandra[1].data
print(chandraData1['PI'])

print("-------------------------------------------FILTERED ENERGIES --------------------------------------------")
chandraData2 = chandraData1['PI'].astype(float)
chandraData2 = chandraData2 / 100
pha_values = chandraData2
x_values = chandraData1['X']
y_values = chandraData1['Y']

# Use numpy.where to filter the data based on energy range criteria
energy_mask = np.where((pha_values >= 1.5) & (pha_values <= 6))
energies = pha_values[energy_mask]
x_positions = x_values[energy_mask]
y_positions = y_values[energy_mask]
result = np.concatenate((energies.reshape(-1,1),x_positions.reshape(-1,1),y_positions.reshape(-1,1)),axis =1)
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

pixel_arcsecond=2.36
maxRadius = 600 #arcsec
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
            radius_intensity[radius]=math.log(radius_intensity[radius]/calculate_area(0,radius))
        else:
            radius_intensity[radius]=math.log(radius_intensity[radius]/calculate_area(radius-5,radius))
#----------------------------------------------------------------------------------------------------------------



#---------------------------------Energy Intensity----------------------------------------------
photonIntensities = {}

#Create 0.2 keV wide intervals
for i in range(30, 121):
    photonIntensities[i/20] = 0

#Round energies to their 1/100 decimals, match them to corresponding energy levels within
for i in result:
    if pixel_arcsecond*(((i[1]-center_point.x)**2+(i[2]-center_point.y)**2)**(0.5)) < maxRadius:
        energy = round(i[0]*20)/20
        for j in photonIntensities.keys():
            if j==energy:
                photonIntensities[j] += 1    

intensityX= []
intensityY = []
for i, j in photonIntensities.items():
    intensityX.append(i)
    intensityY.append(j)
#----------------------------------------------------------------------------------------------



photon_intensity(minRadius, maxRadius, center_point)
x = list(radius_intensity.keys())
y = list(radius_intensity.values())

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))

ax1.plot(intensityX, intensityY, color="darkgreen")
ax2.plot(x, y, color="gold")

ax1.set_title("Photon Count Within " + str(maxRadius) + " Arcsec Radius (1.5, 6 keV)")
ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel('Photon Count')


ax2.set_title('Photon Intensity' + " (" + str(1.5) + ", " + str(6) + " keV)")
ax2.set_xlabel('Distance (Arcsec)')
ax2.set_ylabel('Intensity (log(count/arcsec^2))')


plt.tight_layout()  # Improves spacing between elements
plt.show()
