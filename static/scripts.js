document.addEventListener('DOMContentLoaded', function() {
    // Event listener for the tab buttons
    document.querySelectorAll('.tab button').forEach(function(btn) {
        btn.addEventListener('click', function(e) {
            var tabName = this.getAttribute('data-tab');
            openTab(e, tabName);
        });
    });

    // Function to open a tab
    function openTab(evt, tabName) {
        var i, tabcontent, tablinks;
        // Hide all elements with class="tabcontent" by default
        tabcontent = document.getElementsByClassName('tabcontent');
        for (i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = 'none';
        }

        // Remove the background color of all tablinks/buttons
        tablinks = document.getElementsByClassName('tablinks');
        for (i = 0; i < tablinks.length; i++) {
            tablinks[i].classList.remove('active');
        }

        // Show the specific tab content
        document.getElementById(tabName).style.display = 'block';
        evt.currentTarget.classList.add('active');
    }

    // Simulate a click on the first button/tab of each tab group
    const firstTabButtons = document.querySelectorAll('.tab button:first-child');
    firstTabButtons.forEach(button => button.click());
});
