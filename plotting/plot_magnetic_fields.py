print('PLOTTING PARKER SPIRAL')
theta = (Omega_star / v_sw) * (d_orb - R_star)
x_axis_magnetic_field = d_orb * np.cos(theta)/R_star
y_axis_magnetic_field = d_orb * np.sin(theta)/R_star
plt.figure(figsize=(8, 8))
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax.plot(x_axis_magnetic_field, y_axis_magnetic_field, color='black', lw=2, label="Parker spiral magnetic field line")
circle= plt.Circle((r_orb/R_star, 0), 0.1, color='black', fill=False, linewidth=0.2)
ax.add_patch(circle)
star= plt.Circle((0, 0), 1, color='orange', fill=True, linewidth=2)
ax.add_patch(star)
ax.set_xlabel("X ($R_{star}$)")
ax.set_ylabel("Y ($R_{star}$)")
ax.set_facecolor("white")
ax.text(0, -2,starname,ha='center',fontsize=13)
ax.text(r_orb/R_star, -1,Exoplanet,ha='center',fontsize=13)
ax.axis('equal')

plt.legend()

print(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf')
plt.savefig(FOLDER + '/' + str(Exoplanet.replace(" ", "_")) +'-'+'Parker-spiral-plot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'.pdf', bbox_inches='tight')


