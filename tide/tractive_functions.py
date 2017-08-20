#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 12:44:00 2017

@author: pm7
"""

import numpy as np

def get_xyz_mesh(lon, lat, r):
    # makes a Cartesian mesh on a sphere
    # lon, lat are numpy vectors
    # r is a scalar (the radius)
    x = r * np.outer(np.cos(lon), np.cos(lat))
    y = r * np.outer(np.sin(lon), np.cos(lat))
    z = r * np.outer(np.ones(np.size(lon)), np.sin(lat))
    return x,y,z

def make_arr(a):
    # make sure scalars are numpy arrays
    try:
        len(a)
    except TypeError:
        a=np.array([a])
    return a

def get_xyz(lon, lat, r):
    # makes x,y,z vectors from lon, lat, r
    lon = make_arr(lon)
    lat = make_arr(lat)
    r = make_arr(r)
    x = r * np.cos(lon) * np.cos(lat)
    y = r * np.sin(lon) * np.cos(lat)
    z = r * np.sin(lat)
    return x,y,z
    
def get_R(thx, thy, thz):
    # Cartesian rotation matrix
    c = np.cos(thx)
    s = np.sin(thx)
    Rx = np.matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
    c = np.cos(thy)
    s = np.sin(thy)
    Ry = np.matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    c = np.cos(thz)
    s = np.sin(thz)
    Rz = np.matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return Rz*Ry*Rx
    
def new_xyz(X, Y, Z, R):
    v = np.stack((X, Y, Z), axis=0)
    a = R * v
    # for some reason we have to go through all these
    # shenanigans to get the line or point to plot right
    Xr = make_arr(np.array(a[0,:]).squeeze())
    Yr = make_arr(np.array(a[1,:]).squeeze())
    Zr = make_arr(np.array(a[2,:]).squeeze())
    return Xr, Yr, Zr
    
def get_tractive_scale():
    r = 384000e3 # distance from center of earth to center of moon [m]
    r_e = 6371e3 # radius of the earth [m]
    M = 7.35e22 # mass of the moon [kg]
    E = 5.972e24 # mass of the earth [kg]
    g = 9.8 # acceleration of gravity [m s-2]
    a = g*(M/E)*(r_e/r)**3 # scale of the tractive force [m s-2]
    return a

def get_TF(a, l, d, mlon_rad):
    # returns vectors of the components of the tractive force
    # over the course of a lunar day
    # at latitude l [rad], for lunar declination d [rad]
    # longitude Z (as a way of telling time in lunar days)
    # where Z=0 at the face away from the moon, and Z=pi
    # at the face toward the moon.
    Z = np.pi - mlon_rad
    TFx = (3/2)*a*( np.sin(l)*np.sin(2*d)*np.sin(Z)
                   - np.cos(l)*np.cos(d)**2*np.sin(2*Z) )
    TFy = (3/2)*a*( -(1/2)*np.sin(2*l)*(1 - 3*np.sin(d)**2)
                   - np.cos(2*l)*np.sin(2*d)*np.cos(Z)
                   - (1/2)*np.sin(2*l)*np.cos(d)**2*np.cos(2*Z))
    return TFx, TFy

def draw_vector_ring(ax, mdec_deg, mlon_rad, r, direction):
    if direction=='toward_moon':
        color='g'
    elif direction=='away_from_moon':
        mdec_deg = -mdec_deg
        mlon_rad = mlon_rad + np.pi
        color='gray'
    # ring of vectors facing the Moon
    lon = np.linspace(-np.pi, np.pi, 20)
    lat = np.deg2rad(45) * np.ones_like(lon)
    X, Y, Z = get_xyz(lon, lat, 1.1*r)
    R = get_R(0, np.deg2rad(90-mdec_deg), mlon_rad)
    Xr, Yr, Zr = new_xyz(X, Y, Z, R)
    lon = np.linspace(-np.pi, np.pi, 20)
    lat = np.deg2rad(55) * np.ones_like(lon)
    X, Y, Z = get_xyz(lon, lat, 1.1*r)
    R = get_R(0, np.deg2rad(90-mdec_deg), mlon_rad)
    Xr2, Yr2, Zr2 = new_xyz(X, Y, Z, R)
    ax.quiver(Xr, Yr, Zr, Xr2-Xr, Yr2-Yr, Zr2-Zr,
        linewidth=2, arrow_length_ratio=.5, color=color)
        
def draw_sphere(ax, r):
    lat = np.linspace(-np.pi/2, np.pi/2, 20)
    # first part of the background sphere
    lon = np.linspace(-np.pi, 0, 20)
    x, y, z = get_xyz_mesh(lon, lat, r)
    ax.plot_surface(x, y, z,
                    rstride=1, cstride=1,
                    color='b', linewidth=0, shade=True,
                    alpha=.3)
    # # second part of the background sphere
    lon = np.linspace(0, np.pi, 20)
    x, y, z = get_xyz_mesh(lon, lat, r)
    ax.plot_surface(x, y, z,
                    rstride=1, cstride=1,
                    color='y', linewidth=0, shade=True,
                    alpha=.3)
    # equator line
    lon = np.linspace(-np.pi, np.pi)
    lat = 0*lon
    X, Y, Z = get_xyz(lon, lat, r)
    ax.plot(X, Y, Z, '-g', alpha=.6)
    
def draw_moon(ax, r, mlon_rad, mdec_rad):
    x0, y0, z0 = get_xyz(mlon_rad, mdec_rad, 1.5*r)
    ax.plot(x0, y0, z0, 'og')
    ax.plot([0, x0[0]], [0, y0[0]], [0, z0[0]], '-g')
    #
    lat = np.linspace(-np.pi/2, np.pi/2, 20)
    # first part of the background sphere
    lon = np.linspace(-np.pi, np.pi, 40)
    x, y, z = get_xyz_mesh(lon, lat, r/6)
    ax.plot_surface(x0 + x, y0 + y, z0 + z,
                    rstride=1, cstride=1,
                    color='g', linewidth=0, shade=True,
                    alpha=1)
    
def draw_tractive_pretzels(ax, mdec_rad, mlon_rad_vec, mlon_rad, r):
    # draw the patterns of tractive force over a lunar day
    A = get_tractive_scale()
    AA = A/2
    for tf_lat in np.deg2rad([0, 30, 60]):
        TFx, TFy = get_TF(A, tf_lat, mdec_rad, mlon_rad_vec)
        # beta is an angle that is zero at the north pole
        # and increases to the south (radians)
        beta = np.pi/2 - tf_lat
        dy = (TFx/AA)
        dx = -np.cos(beta) * (TFy/AA)
        dz = np.sin(beta) * (TFy/AA)
        x0, y0, z0 = get_xyz(0, tf_lat, r)
        X = x0 + dx
        Y = y0 + dy
        Z = z0 + dz
        ax.plot(X, Y, Z, '-k')
        # Add a time marker and vector
        TFx, TFy = get_TF(A, tf_lat, mdec_rad, mlon_rad)
        # beta is an angle that is zero at the north pole
        # and increases to the south (radians)
        beta = np.pi/2 - tf_lat
        dy = (TFx/AA)
        dx = -np.cos(beta) * (TFy/AA)
        dz = np.sin(beta) * (TFy/AA)
        X = x0 + dx
        Y = y0 + dy
        Z = z0 + dz
        ax.plot(x0, y0, z0, '*r', markersize=10)
        #ax.plot([x0[0],X], [y0[0],Y], [z0[0],Z], '-r')
        ax.quiver(x0, y0, z0, dx, dy, dz,
            linewidth=2, color='r')

def add_axes(ax, r):
    # labeled coordinate axes
    ax.plot([0, r], [0, 0], [0, 0], '-r')
    ax.text(1.2*r, 0, 0, 'X', color='r', fontweight='bold', fontsize=20)
    ax.plot([0, 0], [0, r], [0, 0], '-b')
    ax.text(0, 1.2*r, 0, 'Y', color='b', fontweight='bold', fontsize=20)
    ax.plot([0, 0], [0, 0], [0, r], '-g')
    ax.text(0, 0, 1.2*r, 'Z', color='g', fontweight='bold', fontsize=20)
    
def make_movie(outdir):
    # and make a movie
    import subprocess
    cmd = ['ffmpeg',
        '-r', '24', # framerate fps
        '-i', outdir + 'plot_%04d.png',
        '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '25',
        outdir + 'movie.mp4']
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if False:
        print('\n-main: screen output from subprocess-')
        print(proc.stdout.decode())
        print('\n-main: errors from subprocess-')
        # for some reason the ffmpeg output ends up in stderr
        print(proc.stderr.decode())
