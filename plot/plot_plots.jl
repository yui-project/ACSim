using Plots	# ここの読み込みにとても時間がかかる
gr()		# バックエンドの設定

# ２要素スカラーのプロット
function plot_2scalar(time,value,filename)
	plot(time,value)
	savefig("./figs/" * filename * ".png")
	
end

function plot_2scalarp(time,value,filename)
	plot(time,value)
	savefig("./figs/" * filename * ".png", marker=:auto)
	
end

function plot_2scalar_range(time,value,xrange, yrange,filename)
	plot(time,value, xlims=(xrange[1], xrange[2]), ylims=(yrange[1], yrange[2]), marker = :auto)
	savefig("./figs/" * filename * ".png")
	
end

# ベクトルを3D上にプロット
function plot_vec(value,filename)
	plot(value[:,1], value[:, 2], value[:, 3])
	xlabel!("x")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/" * filename * ".png")
end

function plot_vecs(value1, value2, valuename, filename)
	plot(value1[:,1], value1[:, 2], value1[:, 3], label = valuename[1])
	plot!(value2[:,1], value2[:, 2], value2[:, 3], label = valuename[2])
	xlabel!("x")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/" * filename * ".png")
end

function plot_vec_range(value, xrange, yrange,zrange, filename)
	plot(value[:,1], value[:, 2], value[:, 3], xlims=(xrange[1], xrange[2]), ylims=(yrange[1], yrange[2]), zlims=(zrange[1], zrange[2]))
	xlabel!("x")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/" * filename * ".png")
end

function plot_2vec_time(time,value,filename)
	plot(time, value[:, 1], value[:, 2])
	xlabel!("time")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/" * filename * ".png")
end

function plot_ecef(x_ecef_log)
	plot(x_ecef_log[:,1], x_ecef_log[:, 2], x_ecef_log[:, 3])
	xlabel!("x")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/" * "x_cecf.png")
end

function plot_3scalar(time, value1, value2, value3, valuename, filename)
	plot(time, value1, label = valuename[1])
	plot!(time, value2, label = valuename[2])
	plot!(time, value3, label = valuename[3])

	savefig("./figs/" * filename * ".png")
end
