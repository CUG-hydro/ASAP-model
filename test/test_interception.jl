

# 测试基本截留计算
@testset "基本截留功能" begin
    minpprate = 0.01
    lai = 3.0
    intercepmax = 0.2 * lai  # 0.6 mm

    # 情况1: 降水量超过截留容量不足量，降水量大于阈值
    precip = 2.0  # mm
    intercepstore = 0.2  # mm
    pet_i = 0.5  # mm

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    deficit = intercepmax - intercepstore  # 0.4 mm
    @test ppdrip ≈ precip - deficit  # 1.6 mm
    @test et_i ≈ 0.0  # 降水量大于阈值，无蒸发损失
    @test new_intercepstore ≈ intercepmax  # 0.6 mm
end

@testset "小降水情况" begin
    lai = 3.0
    intercepmax = 0.2 * lai  # 0.6 mm

    # 情况2: 降水量超过截留容量不足量，但降水量小于阈值
    precip = 0.005  # mm，小于阈值
    intercepstore = 0.2  # mm
    pet_i = 0.5  # mm

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    deficit = intercepmax - intercepstore  # 0.4 mm
    @test ppdrip ≈ 0.0  # -0.395 mm，但逻辑上应该是正值
    @test et_i ≈ 0.205  # 0.5 mm，但受intercepmax限制
    @test new_intercepstore ≈ 0.0
end


@testset "降水不超过不足量" begin
    minpprate = 0.01
    lai = 2.0
    intercepmax = 0.2 * lai  # 0.4 mm

    # 情况3: 降水量不超过截留容量不足量，降水量大于阈值
    precip = 0.1  # mm，大于阈值
    intercepstore = 0.2  # mm
    pet_i = 0.3  # mm

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    @test ppdrip ≈ 0.0
    @test et_i ≈ 0.0  # 降水量大于阈值，无蒸发损失
    @test new_intercepstore ≈ intercepstore + precip  # 0.3 mm
end

@testset "小降水且不超过不足量" begin
    minpprate = 0.01
    lai = 2.0
    intercepmax = 0.2 * lai  # 0.4 mm

    # 情况4: 降水量不超过截留容量不足量，降水量小于阈值
    precip = 0.005  # mm，小于阈值
    intercepstore = 0.2  # mm
    pet_i = 0.1  # mm

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    @test ppdrip ≈ 0.0
    @test et_i ≈ min(intercepstore + precip, pet_i)  # min(0.205, 0.1) = 0.1
    @test new_intercepstore ≈ intercepstore + precip - et_i  # 0.105 mm
end

@testset "边界条件测试" begin
    minpprate = 0.01
    lai = 0.0  # 无叶面积
    intercepstore = 0.0
    pet_i = 1.0
    precip = 5.0

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    @test ppdrip ≈ precip  # 全部穿透
    @test et_i ≈ 0.0
    @test new_intercepstore ≈ 0.0
end

@testset "高叶面积情况" begin
    minpprate = 0.01
    lai = 10.0  # 高叶面积指数
    intercepmax = 0.2 * lai  # 2.0 mm
    intercepstore = 0.5  # mm
    pet_i = 0.8  # mm
    precip = 1.0  # mm

    ppdrip, et_i, new_intercepstore = interception(precip, lai, intercepstore, pet_i)

    deficit = intercepmax - intercepstore  # 1.5 mm
    @test ppdrip ≈ 0.0  # 降水量不超过不足量
    @test et_i ≈ 0.0  # 降水量大于阈值
    @test new_intercepstore ≈ intercepstore + precip  # 1.5 mm
end
