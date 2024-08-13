function loadContent(sectionId) {
  // 隐藏所有内容
  const sections = document.querySelectorAll('.section');
  sections.forEach(section => section.style.display = 'none');

  // 显示指定的内容
  const targetSection = document.getElementById(`${sectionId}-content`);
  if (targetSection) {
      targetSection.style.display = 'block';
  }
}

// 初始化时显示第一个内容
document.addEventListener('DOMContentLoaded', function() {
  loadContent('about');
});




