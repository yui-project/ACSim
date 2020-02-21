using Plots	# ここの読み込みにとても時間がかかる
gr()		# バックエンドの設定

# ２要素スカラーのプロット
function plot_2scalar(time,value,filename)
	plot(time,value)
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
