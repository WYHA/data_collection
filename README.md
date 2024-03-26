**背景**：天文学家通常对一些天文相关的python包有使用需求，比如根据条件进行相关数据检索等等

**目标**：支持用户用语言表述的方式使用python包

**例子**

任务的自然语言表述：查询在Antennae星系周围14角分钟内的相关数据
依赖python包--simbad
我们做的事情：根据用户的任务需求和依赖的python包信息生成python代码-->执行python代码-->返回给用户table数据如下
![](https://cdn.nlark.com/yuque/0/2024/png/29422557/1711414686837-1d3dfb3b-8478-4143-a9c2-a9ff08acae6a.png#averageHue=%23f3f3f3&clientId=u5e323e10-fd54-4&from=paste&id=u122a4e41&originHeight=603&originWidth=1045&originalType=url&ratio=1.5&rotation=0&showTitle=false&status=done&style=none&taskId=ud93c39f3-66fa-4fa9-a895-855e37d1821&title=)

**TODO**
先从几个常用的天文python包入手：astropy，astroquery（如果还有别的常用python包欢迎补充）。希望收集一些任务需求来制作这个方向的数据集，数据收集需求如下：

**数据收集**
需要收集的内容：任务的自然语言表述、依赖的python包、Python code

| 任务的自然语言表述 | 依赖的python包 | Python code | 添加人 |
| --- | --- | --- | --- |
| 用自然语言的方式描述需要astropy，astroquery等完成的任务，组织成一个明确的文本（比如数据检索任务，文本应包含所有的检索条件），如有依赖的文件，可以直接上传，点击左上角插入->文件->本地文件 | 完成任务需要的python包 | 完成任务对应的python code，可以直接将代码复制粘贴至此列。由于太长的python code很难用文本描述清楚，最好将python code的长度控制在100行以内。请将参数正确填充且可以执行成功的python code填写在canonical_solution函数里，将结果通过return进行返回。```def canonical_solution():  ...(此处省略完成任务对应的python code)  return data```| 请务必填写 |
| **例子**：查询在Antennae星系周围14角分钟内的相关数据 | **例子**：astropy | **例子**：```
def canonical_solution():
  import astropy
  ...(此处省略)
  return data
```
 | 王雨菡 |
| 查询glimpse红外巡天星表中以给定txt文件中天体坐标为中心的2角分方形区域内的所有红外天体，txt地址在'/home/usr/object_coordinates.txt'，每一行包含一个天体坐标
**依赖文件**：
 | astropy、astroquery | ```
def canonical_solution():
  from astroquery.ipac.irsa import Irsa
  import astropy.units as u
  Irsa.ROW_LIMIT = 10000
  try:
    with open('/home/usr/object_coordinates.txt', 'r') as file:
      object_coordinates_list = file.readlines()
  except:
      return 'File not found.'
  table = Irsa.query_region(object_coordinates_list, catalog="glimpse_s07", spatial='Box', width=2*u.arcmin).to_pandas()
  return table
```
 | 张天惟 |
| 查询glimpse红外巡天星表中以13:16:43.64 -62:58:31.39坐标为中心的2角分方形区域内的所有红外天体 | astroquery | ```
def canonical_solution():
  from astroquery.ipac.irsa import Irsa
  import astropy.units as u
  Irsa.ROW_LIMIT = 10000
  table = Irsa.query_region("13:16:43.64 -62:58:31.39", catalog="glimpse_s07", spatial='Box', width=2*u.arcmin).to_pandas()
  return table
```
 | 张天惟 |
| 将赤道坐标10.625，41.2转换为银道坐标 | astropy | ```
def canonical_solution():
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')
  return c.galactic
```
 | [https://docs.astropy.org/en/stable/coordinates/](https://docs.astropy.org/en/stable/coordinates/) |
| 将赤道坐标10.625，41.2转换为时分秒写法 | astropy | ```
def canonical_solution():
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')
  c.to_string('hmsdms')
  return c
```
 |  |
| 计算ICRS坐标系下的坐标10，9和FK5坐标系下的坐标11，10之间的空间距离 | astropy | ```
def canonical_solution():
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
  c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='fk5')
  return c1.separation(c2)
```
 |  |
| 将1970年的赤道坐标02h31m49.09s +89d15m50.8s修正为2000年的赤道坐标 | astropy | ```
def canonical_solution():
  from astropy.time import Time
  from astropy.coordinates import SkyCoord
  from astropy.coordinates import FK5
  fk5c = SkyCoord('02h31m49.09s', '+89d15m50.8s',               frame=FK5(equinox=Time('J1970')))
  fk5_2000 = FK5(equinox=Time(2000, format='jyear'))
  return fk5c.transform_to(fk5_2000)
```
 |  |
| 从'tutorials/FITS-images/HorseHead.fits'下载fits文件，做出灰度图和colorbar并保存为horsehead.png | astropy, matplotlib | ```
def canonical_solution():
  import matplotlib.pyplot as plt
  from astropy.visualization import astropy_mpl_style
  from astropy.io import fits
  from astropy.utils.data import get_pkg_data_filename
  plt.style.use(astropy_mpl_style)
  image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
  image_data = fits.getdata(image_file, ext=0)
  plt.figure()
  plt.imshow(image_data, cmap='gray')
  plt.colorbar()
  plt.savefig('horsehead.png')
  return
```
 | [https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html](https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html) |
| 将地址为'HorseHead.fits'的fits文件做出灰度图，横轴纵轴为真实的天体坐标，加网格线，并保存为horsehead.png | astropy, matplotlib | ```
def canonical_solution():
  from astropy.wcs import WCS
  from astropy.io import fits
  import matplotlib.pyplot as plt
  hdu = fits.open('HorseHead.fits')[0]
  wcs = WCS(hdu.header)
  plt.subplot(projection=wcs) 
  plt.imshow(hdu.data, origin='lower', cmap='gray') 
  plt.grid(color='white', ls='solid')
  plt.xlabel('RA (J2000)')
  plt.ylabel('DEC (J2000)')
  plt.savefig('horsehead.png')
  return
```
 |  |
| 将地址为'HorseHead.fits'的fits文件做出灰度图，在最强发射80%的范围内按对数画出5条等高线，横轴纵轴为真实的天体坐标，并保存为horsehead.png | astropy, matplotlib, numpy | ```
# def canonical_solution():
#   from astropy.wcs import WCS
#   from astropy.io import fits
#   import numpy as np
#   import matplotlib.pyplot as plt
#   hdu = fits.open('HorseHead.fits')[0]
#   wcs = WCS(hdu.header)
#   plt.subplot(projection=wcs) 
#   plt.imshow(hdu.data, origin='lower', cmap='gray') 
#   max_v = np.nanmax(hdu.data)
#   min_v = np.nanmin(hdu.data)
#   plt.contour(hdu.data,levels=np.logspace(min_v,max_v*0.8,5),color='red')
#   plt.xlabel('RA (J2000)')
#   plt.ylabel('DEC (J2000)')
#   plt.savefig('horsehead.png')
#   return
def canonical_solution():
  from astropy.wcs import WCS
  from astropy.io import fits
  import numpy as np
  import matplotlib.pyplot as plt
  hdu = fits.open('HorseHead.fits')[0]
  wcs = WCS(hdu.header)
  plt.subplot(projection=wcs) 
  plt.imshow(hdu.data, origin='lower', cmap='gray') 
  max_v = np.nanmax(hdu.data)
  min_v = np.nanmin(hdu.data)
  plt.contour(hdu.data,levels=np.logspace(min_v,max_v*0.8,5))
  plt.xlabel('RA (J2000)')
  plt.ylabel('DEC (J2000)')
  plt.savefig('horsehead.png')
  return
```
 | [https://astropy-astrofrog.readthedocs.io/en/latest/visualization/wcsaxes/images_contours.html](https://astropy-astrofrog.readthedocs.io/en/latest/visualization/wcsaxes/images_contours.html) |
| get the rest frequencies of a list of interstellar molecules? For example, the databases I am planning to use include NIST, splatalogue, or the ALMA database. | astropy, astroquery | ```
def canonical_solution():
  from astropy import units as u
  from astroquery.nist import Nist

  mol_name = 'CO'
    result = Nist.query(mol_name, linename='*', wavelength_range=[0*u.GHz, 1000*u.GHz])

  # Print the results
  for row in result:
      print(row['Wavelength(frequency)'], row['EinsteinA'], row['UpperLevel'], row['LowerLevel'])
```
 | 蒋雪健 |
| concatenate two astropy tables，tables内容分别为{'col1': [1, 2, 3], 'col2': ['a', 'b', 'c']}和{'col1': [4, 5, 6], 'col2': ['d', 'e', 'f']} | astropy | ```
def canonical_solution():
  from astropy.table import Table, vstack

# Create two example tables
  table1 = Table({'col1': [1, 2, 3], 'col2': ['a', 'b', 'c']})
  table2 = Table({'col1': [4, 5, 6], 'col2': ['d', 'e', 'f']})

# Concatenate the tables
  concatenated_table = vstack([table1, table2])

# Print the concatenated table
  print(concatenated_table)
```
 |  |
| 画出fits图像的伪彩色图；叠加另一个fits图像为等高线图；叠加天梯列表为红圈；保存为myfirstplot.png | aplpy | ```
import aplpy
import numpy

gc = aplpy.FITSFigure('fits/2MASS_k.fits')
gc.show_rgb('graphics/2MASS_arcsinh_color.png')

gc.tick_labels.set_font(size='small')

gc.show_contour('fits/mips_24micron.fits', colors='white')

ra, dec = numpy.loadtxt('data/yso_wcs_only.txt', unpack=True)

gc.show_markers(ra, dec, layer='marker_set_1', edgecolor='red',
                facecolor='none', marker='o', s=10, alpha=0.5)

gc.save('myfirstplot.png')
```
 | [Beginner Tutorial — aplpy v2.0.2](https://aplpy.readthedocs.io/en/stable/fitsfigure/quickstart.html) |

可以直接在上表中新增行填写
# data_collection
