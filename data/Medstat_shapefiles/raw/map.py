{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "3dffacdb",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'0.10.2'"
      ]
     },
     "execution_count": 1,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "import geopandas as gpd\n",
    "gpd.__version__"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "b17781b6",
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "bd83dce2",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Requirement already satisfied: folium in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (0.0.0)\n",
      "Requirement already satisfied: requests in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from folium) (2.26.0)\n",
      "Requirement already satisfied: branca>=0.3.0 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from folium) (0.4.2)\n",
      "Requirement already satisfied: numpy in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from folium) (1.21.3)\n",
      "Requirement already satisfied: jinja2>=2.9 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from folium) (3.0.2)\n",
      "Requirement already satisfied: MarkupSafe>=2.0 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from jinja2>=2.9->folium) (2.0.1)\n",
      "Requirement already satisfied: charset-normalizer~=2.0.0 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from requests->folium) (2.0.0)\n",
      "Requirement already satisfied: certifi>=2017.4.17 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from requests->folium) (2021.10.8)\n",
      "Requirement already satisfied: idna<4,>=2.5 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from requests->folium) (3.1)\n",
      "Requirement already satisfied: urllib3<1.27,>=1.21.1 in ./opt/anaconda3/envs/geo_env/lib/python3.9/site-packages (from requests->folium) (1.26.7)\n",
      "Note: you may need to restart the kernel to use updated packages.\n"
     ]
    }
   ],
   "source": [
    "pip install folium"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "b1fc6b5a",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\u001b[34m.\u001b[m\u001b[m/                               .npki_pkcs11.cnf\r\n",
      "\u001b[34m..\u001b[m\u001b[m/                              \u001b[34m.r\u001b[m\u001b[m/\r\n",
      "\u001b[34m.AnySign4PC\u001b[m\u001b[m/                     \u001b[34m.rstudio-desktop\u001b[m\u001b[m/\r\n",
      ".CFUserTextEncoding              \u001b[34m.softforum\u001b[m\u001b[m/\r\n",
      ".DS_Store                        \u001b[34m.ssh\u001b[m\u001b[m/\r\n",
      ".RData                           .tcshrc\r\n",
      ".Rapp.history                    .xonshrc\r\n",
      ".Rhistory                        .zprofile\r\n",
      "\u001b[34m.Trash\u001b[m\u001b[m/                          .zsh_history\r\n",
      ".Xauthority                      \u001b[34m.zsh_sessions\u001b[m\u001b[m/\r\n",
      "\u001b[34m.anaconda\u001b[m\u001b[m/                       .zshrc\r\n",
      ".anyconnect                      7월 공정보고 회의록.hwp\r\n",
      ".bash_profile                    \u001b[34mApplications\u001b[m\u001b[m/\r\n",
      "\u001b[34m.cache\u001b[m\u001b[m/                          \u001b[34mDesktop\u001b[m\u001b[m/\r\n",
      "\u001b[34m.cisco\u001b[m\u001b[m/                          \u001b[34mDocuments\u001b[m\u001b[m/\r\n",
      "\u001b[34m.conda\u001b[m\u001b[m/                          \u001b[34mDownloads\u001b[m\u001b[m/\r\n",
      ".condarc                         Graduate Seminar week 1.txt\r\n",
      "\u001b[34m.config\u001b[m\u001b[m/                         \u001b[34mLibrary\u001b[m\u001b[m/\r\n",
      "\u001b[34m.continuum\u001b[m\u001b[m/                      MEDSTAT_AREAS_2019.shx\r\n",
      "\u001b[34m.cups\u001b[m\u001b[m/                           \u001b[34mMovies\u001b[m\u001b[m/\r\n",
      ".delfino.conf                    \u001b[34mMusic\u001b[m\u001b[m/\r\n",
      "\u001b[34m.fontconfig\u001b[m\u001b[m/                     \u001b[34mPictures\u001b[m\u001b[m/\r\n",
      "\u001b[34m.idlerc\u001b[m\u001b[m/                         \u001b[34mPublic\u001b[m\u001b[m/\r\n",
      "\u001b[34m.ipynb_checkpoints\u001b[m\u001b[m/              \u001b[34mPycharmProjects\u001b[m\u001b[m/\r\n",
      "\u001b[34m.ipython\u001b[m\u001b[m/                        Untitled.ipynb\r\n",
      "\u001b[34m.jupyter\u001b[m\u001b[m/                        \u001b[34mopt\u001b[m\u001b[m/\r\n",
      "\u001b[34m.local\u001b[m\u001b[m/                          stat week 1.txt\r\n",
      "\u001b[34m.matplotlib\u001b[m\u001b[m/\r\n"
     ]
    }
   ],
   "source": [
    "ls -a"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "157b36e2",
   "metadata": {},
   "outputs": [
    {
     "ename": "DriverError",
     "evalue": "MEDSTAT_AREAS_2019.shp: No such file or directory",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mCPLE_OpenFailedError\u001b[0m                      Traceback (most recent call last)",
      "\u001b[0;32mfiona/_shim.pyx\u001b[0m in \u001b[0;36mfiona._shim.gdal_open_vector\u001b[0;34m()\u001b[0m\n",
      "\u001b[0;32mfiona/_err.pyx\u001b[0m in \u001b[0;36mfiona._err.exc_wrap_pointer\u001b[0;34m()\u001b[0m\n",
      "\u001b[0;31mCPLE_OpenFailedError\u001b[0m: MEDSTAT_AREAS_2019.shp: No such file or directory",
      "\nDuring handling of the above exception, another exception occurred:\n",
      "\u001b[0;31mDriverError\u001b[0m                               Traceback (most recent call last)",
      "\u001b[0;32m/var/folders/kc/dxy1rb8j3178f3zmy5kyh9xh0000gn/T/ipykernel_34392/2129797574.py\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[0;32m----> 1\u001b[0;31m \u001b[0mgeo_data\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mgpd\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mread_file\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m'MEDSTAT_AREAS_2019.shp'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mencoding\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;34m'cp949'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mtiles\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;34m'None'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m",
      "\u001b[0;32m~/opt/anaconda3/envs/geo_env/lib/python3.9/site-packages/geopandas/io/file.py\u001b[0m in \u001b[0;36m_read_file\u001b[0;34m(filename, bbox, mask, rows, **kwargs)\u001b[0m\n\u001b[1;32m    199\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    200\u001b[0m     \u001b[0;32mwith\u001b[0m \u001b[0mfiona_env\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 201\u001b[0;31m         \u001b[0;32mwith\u001b[0m \u001b[0mreader\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mpath_or_bytes\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mas\u001b[0m \u001b[0mfeatures\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    202\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    203\u001b[0m             \u001b[0;31m# In a future Fiona release the crs attribute of features will\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/opt/anaconda3/envs/geo_env/lib/python3.9/site-packages/fiona/env.py\u001b[0m in \u001b[0;36mwrapper\u001b[0;34m(*args, **kwargs)\u001b[0m\n\u001b[1;32m    406\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0mwrapper\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m*\u001b[0m\u001b[0margs\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    407\u001b[0m         \u001b[0;32mif\u001b[0m \u001b[0mlocal\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_env\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 408\u001b[0;31m             \u001b[0;32mreturn\u001b[0m \u001b[0mf\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m*\u001b[0m\u001b[0margs\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    409\u001b[0m         \u001b[0;32melse\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    410\u001b[0m             \u001b[0;32mif\u001b[0m \u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0margs\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0;36m0\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mstr\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/opt/anaconda3/envs/geo_env/lib/python3.9/site-packages/fiona/__init__.py\u001b[0m in \u001b[0;36mopen\u001b[0;34m(fp, mode, driver, schema, crs, encoding, layer, vfs, enabled_drivers, crs_wkt, **kwargs)\u001b[0m\n\u001b[1;32m    254\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    255\u001b[0m         \u001b[0;32mif\u001b[0m \u001b[0mmode\u001b[0m \u001b[0;32min\u001b[0m \u001b[0;34m(\u001b[0m\u001b[0;34m'a'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'r'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 256\u001b[0;31m             c = Collection(path, mode, driver=driver, encoding=encoding,\n\u001b[0m\u001b[1;32m    257\u001b[0m                            layer=layer, enabled_drivers=enabled_drivers, **kwargs)\n\u001b[1;32m    258\u001b[0m         \u001b[0;32melif\u001b[0m \u001b[0mmode\u001b[0m \u001b[0;34m==\u001b[0m \u001b[0;34m'w'\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/opt/anaconda3/envs/geo_env/lib/python3.9/site-packages/fiona/collection.py\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, path, mode, driver, schema, crs, encoding, layer, vsi, archive, enabled_drivers, crs_wkt, ignore_fields, ignore_geometry, **kwargs)\u001b[0m\n\u001b[1;32m    160\u001b[0m             \u001b[0;32mif\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mmode\u001b[0m \u001b[0;34m==\u001b[0m \u001b[0;34m'r'\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    161\u001b[0m                 \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msession\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mSession\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 162\u001b[0;31m                 \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msession\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mstart\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    163\u001b[0m             \u001b[0;32melif\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mmode\u001b[0m \u001b[0;32min\u001b[0m \u001b[0;34m(\u001b[0m\u001b[0;34m'a'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'w'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    164\u001b[0m                 \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msession\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mWritingSession\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32mfiona/ogrext.pyx\u001b[0m in \u001b[0;36mfiona.ogrext.Session.start\u001b[0;34m()\u001b[0m\n",
      "\u001b[0;32mfiona/_shim.pyx\u001b[0m in \u001b[0;36mfiona._shim.gdal_open_vector\u001b[0;34m()\u001b[0m\n",
      "\u001b[0;31mDriverError\u001b[0m: MEDSTAT_AREAS_2019.shp: No such file or directory"
     ]
    }
   ],
   "source": [
    "geo_data = gpd.read_file('MEDSTAT_AREAS_2019.shp', encoding='cp949', tiles='None')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b298cffa",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
