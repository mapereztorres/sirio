#print('PLOTTING PARKER SPIRAL')
theta = (Omega_star / v_sw) * (d_orb - R_star)
x = d_orb * np.cos(theta)/R_star
y = d_orb * np.sin(theta)/R_star
plt.figure(figsize=(8, 8))
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax.plot(x, y, color='black', lw=2, label="Parker spiral magnetic field line")


ax.set_xlabel("X ($R_{star}$)")
ax.set_ylabel("Y ($R_{star}$)")
ax.set_facecolor("white")

print(ax.get_ylim())
limy=ax.get_ylim()
limy=limy[0]



ax.set_ylim([-3,3])

plt.legend()

print(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf')
plt.savefig(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf', bbox_inches='tight')


