#r = (v / omega) * theta + r0
#R0=R_STAR
#v=v_sw
#r=d_orb
print('PLOTTING PARKER SPIRAL')
theta = (Omega_star / v_sw) * (d_orb - R_star)
x = d_orb * np.cos(theta)/R_star
y = d_orb * np.sin(theta)/R_star
plt.figure(figsize=(8, 8))
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax.plot(x, y, color='black', lw=2, label="Parker spiral magnetic field line")
#ax.scatter(0, 0, color='yellow', s=(10)**2,label=starname) 
#ax.scatter(r_orb/R_star, 0 , color='black', s=10, label=Exoplanet)  
circle= plt.Circle((r_orb/R_star, 0), 0.1, color='black', fill=False, linewidth=0.2)
ax.add_patch(circle)
star= plt.Circle((0, 0), 1, color='orange', fill=True, linewidth=2)
ax.add_patch(star)
ax.set_xlabel("X ($R_{star}$)")
ax.set_ylabel("Y ($R_{star}$)")
ax.set_facecolor("white")
ax.text(0, -2,starname,ha='center',fontsize=13)
ax.text(r_orb/R_star, -1,Exoplanet,ha='center',fontsize=13)
#plt.title("Espiral de Parker (Vista desde el eje Z)")
#plt.title("Parker spiral")
#plt.legend()
#plt.grid()
#plt.axis("equal")  # Mantiene la escala de los ejes iguales
#plt.savefig('OUTPUT/foo.pdf')
ax.axis('equal')
#ax.set_aspect('equal')
#black_patch = mpatches.Patch(color='black', label='$R_{mp}$')
#red_patch = mpatches.Patch(color='red', label='$R_{eff}$')
#ax.legend(handles=[black_patch,red_patch],loc='upper left',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)

#circulo_legenda = mpatches.Circle((0, 0), radius=0.2, color='blue', fill=False, linewidth=2)
#plt.legend(handles=[circulo_legenda], labels=["Círculo de referencia"], loc="upper right")

plt.legend()

print(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf')
plt.savefig(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf', bbox_inches='tight')
#plt.show()

