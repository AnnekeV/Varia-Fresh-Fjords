import xarray as xr
import matplotlib.pyplot as plt
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

print("opening dataset")
start_time = time.time()


def print_time_elapsed(start_time=start_time):
    elapsed_time = time.time() - start_time
    print(f"--- {elapsed_time:.0f} seconds ---")


base_folder_terminal = "/science/projects/imau02/rapid/Anneke/CARRA/"
base_folder_local = "/Volumes/imau02/rapid/Anneke/CARRA/"
base_folder = base_folder_terminal

ds_tp = xr.open_mfdataset(
    base_folder
    + "Monthly/total_precipitation.CARRA.west_domain.forecast.Monthly.2009-2023.nc"
)
print_time_elapsed()

ds_tp = ds_tp.sortby('time')

print("making the plot")
fig = (
    ds_tp["tp"]
    .squeeze()[:, ::10, ::10]
    .plot(
        x="x",
        y="y",
        col="time",
        col_wrap=12,
        robust=True,
        cmap="viridis",
        figsize=(40, 50),
    )
)
print_time_elapsed()

print("saving the plot")
plt.savefig("plot_quicklook_tp.png")
print_time_elapsed()

print("showing the plot")
plt.show()
